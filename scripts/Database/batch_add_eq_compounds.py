import pickle
import queue
import re
import signal
import time
from itertools import islice
from multiprocessing import Process, Queue, set_start_method
from pathlib import Path

import pandas as pd
import sqlalchemy
from minedatabase.pickaxe import Pickaxe

from typing import List, Tuple, Generator

import sys

# Set up some equilibrator items
# Import and define working database URI
cwd = Path(__file__).parent
URI_EQ = open(cwd / "eq_uris.txt", "r").read().strip("\n")

# Define some utility functions
class TimeOutException(Exception):
    pass

def alarm_handler(signum, frame):
    raise TimeOutException

class Timer:
    """Timer class to keep track of timings.

    Prints out start/stop messages and timing.
    """
    def __init__(self, message: str) -> None:
        self.message = message
        self.start_time = time.time()
        self.checkpoint_time = time.time()

    def start_timer(self, print_message: bool = True):
        self.start_time = time.time()
        if print_message:
            print(f"Starting {self.message}, setting timer to 0")

    def print_time_from_start(self, print_message: bool = True, header: str = "", footer: str = ""):
        if print_message:
            if header:
                print(header)
            print(f"Time for {self.message} from start: {time.time() - self.start_time}")
            if footer:
                print(footer)
        else:
            print(time.time() - self.start_time)

    def set_checkpount(self, print_message: bool = True):
        self.checkpoint_time = time.time()

    def print_time_from_checkpoint(self, print_message: bool = True, header: str = "", footer: str = "", set_new_checkpoint: bool = True):
        if print_message:
            if header:
                print(header)
            print(f"Time for {self.message} from last checkpoint: {time.time() - self.checkpoint_time}")
            if footer:
                print(footer)

        if set_new_checkpoint:
            self.checkpoint_time = time.time()

def yield_chunks_from_pickaxe(pk: Pickaxe, chunk_size: int) -> Generator[List[Tuple[str, str]], None, None]:
    """Provides a generator yielding lists of compound SMILES

        Returns a generator that yiels a list of compound SMILES from a pickaxe object.

        Parameters
        ----------
        pk : Pickaxe
            A pickaxe object, whose compound SMILES will be returned
        chunk_size : int
            The size of the chunks to send to the child processes

        Yields
        ------
        Generator[List[Tuple[str, str]], None, None]
            A generator containing lists of compound SMILES and compound IDs
        """
    #TODO Only using 100 cpds, delete that for more
    white_list = ['C98634908f9c25af583dabc72a1d5a2e26b154f77', 'X85753d94c72a1ae851d34ac6ac8b0c9828cfc734', 'C8f24b2071926f02213b99cb6f7df5e035b9ce890', 'Cec77ea281f69ca989bbef0a4c7794128a196c716', 'C53ebd0b2e0c6e4c57239f8730243822487e19930', 'C8917b48859fb749bc90b2140befd2422c11bbba0', 'C38fa6d3508da96f383e5b60ec41001692bc645cf', 'Xc9b0df14e690ee411edec0aad2286b732e431b96', 'C4d0068f7803b26ab3c9d318cf8abdd4d4a584534', 'C484f65e0d1a503c751d32a31cc561c8d7d786824', 'C50cfd439c30042ba827bc8c0353b45f6ae72f270', 'C6ec1611229ff4fc7a19244967c7716266fc021a1']
    compound_list = [(cpd["SMILES"], cpd["_id"]) for cpd in pk.compounds.values() if cpd["_id"] in white_list]

    for i in range(0, len(compound_list), chunk_size):
        yield compound_list[i : i + chunk_size]

# This is the heavylifter -- Function that is used to initialize children processes
def get_compounds_in_queue(
    compounds_in: queue,
    good_compounds: queue,
    bad_compounds: queue,
    pid: int=0,
    timeout: float=360,
    bypass_chemaxon: bool=True,
    save_empty_compounds: bool=True
) -> None:
    """Add compounds to eQuilibrator from a queue

    Spawns a child process that receives information from a parent. Tries to add them to eQuilibrator and returns bad compounds.

    Parameters
    ----------
    compounds_in : queue
        Queue of input compounds. In the format of [(SMILES, ID), ...]
    bad_compounds : queue
        Queue of compounds that failed. In the format of [(SMILES, ID), ...]
    pid : int, optional
        The process ID given by the parent process., by default 0
    timeout : float, optional
        The timeout for eQuilibrator getting compounds., by default 360
    bypass_chemaxon : bool, optional
        Whether or not to bypass chemaxon if no pKas can be computed, by default True
    save_empty_compounds : bool, optional
        Whether or not to save an empty compound if compound can't be decomposed, by default True
    """
    signal.signal(signal.SIGALRM, alarm_handler)
    process_timer = Timer(f"Child Process {pid} Timer")
    process_timer.start_timer(print_message=True)
    loop = 0
    
    # with open(f"./stdouts/Error_{pid}.txt", "w") as sys.stdout:
    # Load and init
    print("Importing eQuilibrator")
    from equilibrator_assets.local_compound_cache import LocalCompoundCache
    from equilibrator_cache.compound_cache import CompoundCache

    print(f"Connecting {pid} to postgresql DB and creating Local Compound Cache.\n")
    lcp = LocalCompoundCache()
    lcp.ccache = CompoundCache(sqlalchemy.create_engine(URI_EQ))
    process_timer.print_time_from_start(
        print_message=True,
        header=f"{'-'*15}\nChild Process {pid} succesfully initialized",
        footer="-"*15
    )

    # Loop and get compounds
    while True:
        process_timer.set_checkpount(False)
        loop += 1
        # Try to get the list of SMILES
        try:
            compound_list = compounds_in.get(timeout=2)
            process_timer.print_time_from_checkpoint(
                print_message=True,
                header=f"{'-'*15}\nChild Process {pid} | Loop {loop}: Sucessfuly got new compound list.",
                footer="-"*15,
                set_new_checkpoint=False
            )

        except queue.Empty:
            print(
                print(f"Child Process {pid} | Loop {loop}: Found Empty Queue after ")
            )
            process_timer.print_time_from_start(
                print_message=True,
                header=f"{'-'*15}\nChild Process {pid} | Loop {loop}: Empty queue found, terminating.",
                footer="-"*15,
            )
            break

        else: # We now have a list of compounds, grab it from the database
            # Reset timeout timer
            signal.alarm(timeout)
            process_timer.print_time_from_checkpoint()
            print(f"Child Process {pid} | Loop {loop}: Evaluating reactions with timeout {timeout}")
            process_timer.set_checkpount(False)
            
            # Try to get the compounds from the DB
            try:
                # We will filter out compounds with a "*" in them. #TODO what to do about this?
                # We also define our coco accession with pk_{compound_id} for easy lookup later
                compound_list = [(compound_smiles, f"pk_{compound_id}", compound_id) for compound_smiles, compound_id in compound_list if "*" not in compound_smiles]

                # Create a dataframe for equilibrator to add from
                df_compounds_to_add = pd.DataFrame(
                    compound_list, columns=["struct", "coco_id", "name"]
                )

                # Add compounds to equilibrator
                lcp.add_compounds(
                    df_compounds_to_add, save_empty_compounds=save_empty_compounds, bypass_chemaxon=bypass_chemaxon
                )

                process_timer.print_time_from_checkpoint(
                    print_message=True,
                    header=f"Child Process {pid} | Loop {loop}: Succesfully added compounds.",
                    footer=""
                )
                
                good_compounds.put((df_compounds_to_add.struct, df_compounds_to_add.name))

            except TimeOutException:
                process_timer.print_time_from_checkpoint(
                    print_message=True,
                    header=f"Child Process {pid} | Loop {loop}: Timed out adding compounds.",
                    footer=""
                )
                #TODO how do we handle bad compounds? For now just send them back
                bad_compounds.put(compound_list)

            except BaseException as e:
                process_timer.print_time_from_checkpoint(
                    print_message=True,
                    header=f"Child Process {pid} | Loop {loop}: Error: {e}",
                    footer=""
                )
                                    #TODO how do we handle bad compounds? For now just send them back
                bad_compounds.put(compound_list)

            # Turn timeout alarm off until next loop
            signal.alarm(0)

        # Shut everything down and return
    lcp.ccache.session.close()
    del lcp
    process_timer.print_time_from_start(
        print_message=True,
        header=f"{'-'*15}\nChild Process {pid} | Loop {loop}: Succesfully terminated.",
        footer="-"*15
    )
    return True

# This is the parent function that spawns children
def insert_chunk_queue(pk, n_processes, chunk_size, bypass_chemaxon, save_empty_compounds):
    """Insert chunked compounds in queue"""
    children = []
    timeout = 20
    n_compounds = len(pk.compounds)

    # Make queues for communication between parent and child
    compound_list_queue = Queue()
    good_compounds_queue = Queue()
    bad_compounds_queue = Queue()

    # Populate the input queue
    for compound_list in yield_chunks_from_pickaxe(pk, chunk_size):
        compound_list_queue.put(compound_list)

    # Start timer
    parent_timer = Timer("Parent Process")

    # Start processes and set up queues, etc.
    print("Starting Processes")
    for i in range(n_processes):
        p = Process(
            target=get_compounds_in_queue,
            args=(compound_list_queue, good_compounds_queue, bad_compounds_queue, i, timeout, bypass_chemaxon, save_empty_compounds),
        )
        children.append(p)

    for child in children:
        child.start()

    parent_timer.print_time_from_start(
        print_message=True,
        header=f"{'-'*15}\nParent process succesfully spawned children.",
        footer="-"*15
    )
    parent_timer.set_checkpount(False)

    # Specify some other parameters
    # Time to wait to print parent status
    print_dt = 20
    good_compounds = []
    bad_compounds = []

    # Track how many succesful batches there were
    good_i = 0
    bad_i = 0

    # Join + start processes and get outputs
    print("Joining and Starting")
    while any([child.is_alive() for child in children]):
        for child in children:
            # Not sure what this does... but it's here. #TODO figure out what it does
            child.join(timeout=0.1)

        try:
            good_compounds.extend(good_compounds_queue.get_nowait())
            good_i += chunk_size
        except queue.Empty:
            pass

        try:
            bad_compounds.append(bad_compounds_queue.get_nowait())
            bad_i += chunk_size
        except queue.Empty:
            pass

        if time.time() - parent_timer.checkpoint_time >= print_dt:
            parent_timer.print_time_from_checkpoint(
                True,
                header=f"{'-'*15}\nParent Process Checkpoint Update.\n\tProgress: {(good_i + bad_i)} of {n_compounds//chunk_size + (1 if n_compounds%chunk_size else 0)}\n\tGood: {good_i*chunk_size}\n\tBad: {bad_i*chunk_size}",
                footer="-"*15
            )
    # TODO
    # Check how children were closed
    for child in children:
        p.join()

    parent_timer.print_time_from_start(
        True,
        header=f"{'-'*15}\nParent Process Completed.\n\tGood: {good_i*chunk_size}\n\tBad: {bad_i*chunk_size}",
        footer="-"*15
    )

    return good_compounds, bad_compounds

def main():
    #TODO argparse stuff for inputs
    n_processes = 1
    chunk_size = 1
    bypass_chemaxon = True
    save_empty_compounds = True

    print("Loading Pickaxe and Generating eQuilibrator Compounds")

    pk = Pickaxe()
    # pk_location = cwd / "aGDM_example.pk"
    pk_location = "/home/stef/bottle/data/raw_expansions/succinate_to_mvacid_gen_4_tan_sample_1_n_samples_1000.pk"
    pk.load_pickled_pickaxe(pk_location)

    print(len(pk.compounds))
    print("\nLoaded Pickaxe, beginning eQuilibrator compound generation.")

    print("Beginning get Reactions")
    good_compounds, bad_compounds = insert_chunk_queue(pk, n_processes, chunk_size, bypass_chemaxon, save_empty_compounds)

    print("Saving good and bad compounds")
    with open("good_eq_compounds.csv", "w") as f:
        f.writelines([f"{compound_SMILES},{compound_id}" for compound_SMILES, compound_id in good_compounds])
    
    with open("bad_eq_compounds.csv", "w") as f:
        f.writelines([f"{compound_SMILES},{compound_id}" for compound_SMILES, compound_id in bad_compounds])

if __name__ == "__main__":
    # This is important for MacOS/Linux. Make it so each child process is NEW
    # If you don't do this it shares resources, namely the local compound cache
    # This really screws with adding things, as there is now only one compound cache
    set_start_method("spawn")

    # If you don't do spawn the local cache will get shared (potentially, may not actually happen) between the processes
    # This is just a way to keep the children independent

    main()
