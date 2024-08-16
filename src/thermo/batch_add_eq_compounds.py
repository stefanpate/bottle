import os
import pickle
import queue
import re
import signal
import sys
import time
from itertools import islice
from multiprocessing import Process, Queue, set_start_method
from pathlib import Path
from typing import Generator, List, Tuple

import pandas as pd
import sqlalchemy
from minedatabase.pickaxe import Pickaxe

# Set up some equilibrator items
# Import and define working database URI
# cwd = Path(__file__).parent
LOG_PATH = "../logs/batch_add_logs"
URI_EQ = open("../artifacts/eq_uris.uri", "r").read().strip("\n")

def print_flush(*args, flush=True, **kwargs):
    print(*args, **kwargs, flush=flush)

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
            print_flush(f"Starting {self.message}, setting timer to 0")

    def print_time_from_start(
        self, print_message: bool = True, header: str = "", footer: str = ""
    ):
        if print_message:
            if header:
                print_flush(header)
            print_flush(
                f"Time for {self.message} from start: {time.time() - self.start_time}"
            )
            if footer:
                print_flush(footer)
        else:
            print_flush(time.time() - self.start_time)

    def set_checkpount(self, print_message: bool = True):
        self.checkpoint_time = time.time()

    def print_time_from_checkpoint(
        self,
        print_message: bool = True,
        header: str = "",
        footer: str = "",
        set_new_checkpoint: bool = True,
    ):
        if print_message:
            if header:
                print_flush(header)
            print_flush(
                f"Time for {self.message} from last checkpoint: {time.time() - self.checkpoint_time}"
            )
            if footer:
                print_flush(footer)

        if set_new_checkpoint:
            self.checkpoint_time = time.time()


def yield_chunks_from_pickaxe(
    pk: Pickaxe,
    chunk_size: int,
    return_unadded_only: bool = False,
    ignore_asterisk: bool = True,
    max_entries: int = None,
) -> Generator[List[Tuple[str, str]], None, None]:
    """
    Yield chunks of compounds from a Pickaxe object.

    Args:
        pk (Pickaxe): The Pickaxe object containing the compounds.
        chunk_size (int): The size of each chunk.
        return_unadded_only (bool, optional): Whether to return only the unadded compounds.
        max_entries (int, optional): The maximum number of entries to process.

    Yields:
        Generator[List[Tuple[str, str]], None, None]: A generator that yields chunks of compounds.

    """
    compound_list = [
        (cpd["SMILES"], cpd["_id"]) for cpd in pk.compounds.values() if cpd["_id"]
    ]

    if ignore_asterisk:
        compound_list = [
            (compound_SMILES, compound_id)
            for compound_SMILES, compound_id in compound_list
            if "*" not in compound_SMILES
        ]

    if return_unadded_only:
        from src.thermo.local_compound_cache import LocalCompoundCache
        from equilibrator_cache.compound_cache import CompoundCache
        from src.thermo.pickaxe_thermodynamics import PickaxeThermodynamics
        print_flush(
            f"Connecting Parent to postgresql DB and creating Local Compound Cache to see which compounds already exist.\n"
        )
        lcp = LocalCompoundCache()
        lcp.ccache = CompoundCache(sqlalchemy.create_engine(URI_EQ))

        PT = PickaxeThermodynamics(lcp)

        PT.generate_eQ_compound_dict_from_pickaxe(pk)
        already_added  = [
            (compound_SMILES, compound_id)
            for compound_SMILES, compound_id in compound_list
            if compound_id in PT.eQ_compound_dict
        ]

        compound_list = [
            (compound_SMILES, compound_id)
            for compound_SMILES, compound_id in compound_list
            if compound_id not in PT.eQ_compound_dict
        ]

        print_flush(f"Found {len(already_added)} compounds already added to eQuilibrator.")

    if max_entries:
        compound_list = compound_list[:max_entries]

    for i in range(0, len(compound_list), chunk_size):
        yield compound_list[i : i + chunk_size]


# This is the heavylifter -- Function that is used to initialize children processes
def get_compounds_in_queue(
    compounds_in: queue,
    good_compounds: queue,
    bad_compounds: queue,
    pid: int = 0,
    timeout: float = 360,
    bypass_chemaxon: bool = True,
    save_empty_compounds: bool = True,
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
    sys.stdout = open(os.path.join(LOG_PATH, f"stdout_{pid}.log"), "w")
    sys.stderr = open(os.path.join(LOG_PATH, f"stderr_{pid}.log"), "w")

    signal.signal(signal.SIGALRM, alarm_handler)
    process_timer = Timer(f"Child Process {pid} Timer")
    process_timer.start_timer(print_message=True)
    loop = 0

    # with open(f"./stdouts/Error_{pid}.txt", "w") as sys.stdout:
    # Load and init
    print_flush("Importing eQuilibrator")
    from src.thermo.local_compound_cache import \
        LocalCompoundCache
    from equilibrator_cache.compound_cache import CompoundCache

    print_flush(
        f"Connecting {pid} to postgresql DB and creating Local Compound Cache.\n"
    )
    lcp = LocalCompoundCache()
    lcp.ccache = CompoundCache(sqlalchemy.create_engine(URI_EQ))

    print_flush(lcp)
    
    process_timer.print_time_from_start(
        print_message=True,
        header=f"{'-'*15}\nChild Process {pid} succesfully initialized",
        footer="-" * 15,
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
                footer="-" * 15,
                set_new_checkpoint=False,
            )

        except queue.Empty:
            print_flush(
                print_flush(
                    f"Child Process {pid} | Loop {loop}: Found Empty Queue after "
                )
            )
            process_timer.print_time_from_start(
                print_message=True,
                header=f"{'-'*15}\nChild Process {pid} | Loop {loop}: Empty queue found, terminating.",
                footer="-" * 15,
            )
            break

        else:  # We now have a list of compounds, grab it from the database
            # Reset timeout timer
            signal.alarm(timeout)
            process_timer.print_time_from_checkpoint()
            print_flush(
                f"Child Process {pid} | Loop {loop}: Trying to evaluate {len(compound_list)} compounds with timeout {timeout}"
            )
            process_timer.set_checkpount(False)

            # Try to get the compounds from the DB
            try:
                # We also define our coco accession with pk_{compound_id} for easy lookup later
                compound_list = [
                    (compound_smiles, f"pk_{compound_id}", compound_id)
                    for compound_smiles, compound_id in compound_list
                ]
                # Create a dataframe for equilibrator to add from
                df_compounds_to_add = pd.DataFrame(
                    compound_list, columns=["struct", "coco_id", "name"]
                )

                # Add compounds to equilibrator
                # Remove the pickaxe to use regular code
                lcp.add_pickaxe_compounds(
                    df_compounds_to_add,
                    save_empty_compounds=save_empty_compounds,
                    bypass_chemaxon=bypass_chemaxon,
                )

                process_timer.print_time_from_checkpoint(
                    print_message=True,
                    header=f"Child Process {pid} | Loop {loop}: Succesfully added compounds.",
                    footer="",
                )

                good_compounds.put(
                    [*zip(df_compounds_to_add.struct, df_compounds_to_add.name)]
                )



            except TimeOutException:
                process_timer.print_time_from_checkpoint(
                    print_message=True,
                    header=f"Child Process {pid} | Loop {loop}: Timed out adding compounds.",
                    footer="",
                )
                # TODO how do we handle bad compounds? For now just send them back
                bad_compounds.put(compound_list)

            except BaseException as e:
                process_timer.print_time_from_checkpoint(
                    print_message=True,
                    header=f"Child Process {pid} | Loop {loop}: Error: {e}",
                    footer="",
                )
                # Our database is potentially bad... nuke and remake it
                lcp.ccache.session.close()
                lcp = LocalCompoundCache()
                lcp.ccache = CompoundCache(sqlalchemy.create_engine(URI_EQ))
                # TODO how do we handle bad compounds? For now just send them back
                bad_compounds.put(compound_list)

            # Turn timeout alarm off until next loop
            signal.alarm(0)

            # Delete lcp and reconnect
            # process_timer.print_time_from_checkpoint(
            #     print_message=True,
            #     header=f"Child Process {pid} | Loop {loop}: Closing and reconnecting to database.",
            #     footer="",
            # )
            # lcp.ccache.session.close()
            # lcp.ccache = CompoundCache(sqlalchemy.create_engine(URI_EQ))
            # process_timer.print_time_from_checkpoint(
            #     print_message=True,
            #     header=f"Child Process {pid} | Loop {loop}: Succesfully closed and reconnected to database.",
            #     footer="",
            # )

        # Shut everything down and return
    lcp.ccache.session.close()
    del lcp
    process_timer.print_time_from_start(
        print_message=True,
        header=f"{'-'*15}\nChild Process {pid} | Loop {loop}: Succesfully terminated.",
        footer="-" * 15,
    )
    return True


# This is the parent function that spawns children
def insert_chunk_queue(
    pk, n_processes, chunk_size, bypass_chemaxon, save_empty_compounds, return_unadded_only=True, ignore_asterisk=True, max_entries=None
):
    """Insert chunked compounds in queue"""
    children = []
    timeout = 60 * 100  # 4 minutes
    n_compounds = len(pk.compounds)

    # Make queues for communication between parent and child
    compound_list_queue = Queue()
    good_compounds_queue = Queue()
    bad_compounds_queue = Queue()
    

    # Populate the input queue
    n_compounds = 0
    n_chunks = 0
    print_flush("Populating Input Queue.")
    for compound_list in yield_chunks_from_pickaxe(pk=pk, chunk_size=chunk_size, return_unadded_only=return_unadded_only, ignore_asterisk=ignore_asterisk, max_entries=max_entries):
        compound_list_queue.put(compound_list)
        n_compounds += len(compound_list)
        n_chunks += 1
    print_flush(
        f"\tPopulated Input Queue with {n_compounds} compounds in {n_chunks} chunks."
    )

    # Start timer
    parent_timer = Timer("Parent Process")

    # Start processes and set up queues, etc.
    print_flush("Starting Processes from Parent.")
    for i in range(n_processes):
        print(f"\tStarting Process, {i}")
        p = Process(
            target=get_compounds_in_queue,
            args=(
                compound_list_queue,
                good_compounds_queue,
                bad_compounds_queue,
                i,
                timeout,
                bypass_chemaxon,
                save_empty_compounds,
            ),
        )
        children.append(p)

    for child in children:
        child.start()

    parent_timer.print_time_from_start(
        print_message=True,
        header=f"{'-'*15}\nParent process succesfully spawned children.",
        footer="-" * 15,
    )
    parent_timer.set_checkpount(False)

    # Specify some other parameters
    # Time to wait to print parent status
    print_dt = 30
    good_compounds = []
    bad_compounds = []

    # Track how many succesful batches there were
    good_n_chunks = 0
    bad_n_chunks = 0

    good_n_compounds = 0
    bad_n_compounds = 0

    # Join + start processes and get outputs
    print_flush("Joining and Starting")
    while any([child.is_alive() for child in children]):
        for child in children:
            # Not sure what this does... but it's here. #TODO figure out what it does
            child.join(timeout=0.1)

        try:
            good_compounds.extend(good_compounds_queue.get_nowait())
            good_n_chunks += 1
            good_n_compounds = len(good_compounds)

        except queue.Empty:
            pass

        try:
            bad_compounds.extend(bad_compounds_queue.get_nowait())
            bad_n_chunks += 1
            bad_n_compounds = len(bad_compounds)

        except queue.Empty:
            pass

        if time.time() - parent_timer.checkpoint_time >= print_dt:
            parent_timer.print_time_from_checkpoint(
                True,
                header=f"{'-'*15}\nParent Process Checkpoint Update.\n\tProgress:\n\t\tChunks: {(good_n_chunks + bad_n_chunks)} of {n_chunks}\n\t\tCompounds Total: {(good_n_compounds + bad_n_compounds)} of {n_compounds}\n\t\t\tGood: {good_n_compounds}\n\t\t\tBad: {bad_n_compounds}",
            )
            parent_timer.print_time_from_start(True, footer="-" * 15)
    # TODO
    # Check how children were closed
    for child in children:
        p.join()

    parent_timer.print_time_from_start(
        True,
        header=f"{'-'*15}\nParent Process Completed.\n\tGood: {good_n_compounds}\n\tBad: {bad_n_compounds}",
        footer="-" * 15,
    )

    for child in children:
        print_flush(f"Child is alive: {child.is_alive()}")

    return good_compounds, bad_compounds


def add_compounds_to_eQ(pk:Pickaxe):
    # Change these parameters to change how the script runs
    n_processes = 1
    chunk_size = 50
    bypass_chemaxon = True
    save_empty_compounds = True
    return_unadded_only = True
    ignore_asterisk = True
    max_entries = None

     # Make log folder
    os.makedirs(LOG_PATH, exist_ok=True)
    # Remove all files from batch_add_logs folder
    for file in os.listdir(LOG_PATH):
        os.remove(os.path.join(LOG_PATH, file))

    # Open the files
    stdout_file = open(os.path.join(LOG_PATH, "stdout_parent.log"), "w")
    stderr_file = open(os.path.join(LOG_PATH, "stderr_parent.log"), "w")

    # Redirect stdout and stderr
    stdout = sys.stdout
    stderr = sys.stderr
    sys.stdout = stdout_file
    sys.stderr = stderr_file
    
    print_flush("Loading Pickaxe and Generating eQuilibrator Compounds")
    print(f"n_processes: {n_processes}\nchunk_size: {chunk_size}\nreturn_unadded_only: {return_unadded_only}\nignore_asterisk: {ignore_asterisk}")
    print(f"max_entries: {max_entries}")

    print_flush(len(pk.compounds))
    print_flush("\nLoaded Pickaxe, beginning eQuilibrator compound generation.")

    print_flush("Beginning get Reactions")
    good_compounds, bad_compounds = insert_chunk_queue(
        pk, n_processes, chunk_size, bypass_chemaxon, save_empty_compounds, return_unadded_only, ignore_asterisk, max_entries
    )

    print_flush("Saving Bad compounds")
    with open(os.path.join(LOG_PATH, "good_eq_compounds.csv"), "w") as f:
        f.writelines([f"{compound_SMILES},{compound_id}" for compound_SMILES, compound_id in good_compounds])

    with open(os.path.join(LOG_PATH, "bad_eq_compounds.csv"), "w") as f:
        f.writelines([f"{compound_SMILES},{compound_id}" for compound_SMILES, _, compound_id in bad_compounds])

    # Close the files
    stdout_file.close()
    stderr_file.close()

    # Reset stdout stderr
    sys.stdout = stdout
    sys.stderr = stderr


if __name__ == "__main__":
    # This is important for MacOS/Linux. Make it so each child process is NEW
    # If you don't do this it shares resources, namely the local compound cache
    # This really screws with adding things, as there is now only one compound cache
    set_start_method("spawn")
    pk_path = ''
    pk = Pickaxe()
    pk.load_pickled_pickaxe(pk_path)

    # Your code here
    add_compounds_to_eQ(pk)

