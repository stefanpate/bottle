import pandas as pd
from src.config import filepaths
from minedatabase.pickaxe import Pickaxe
from itertools import product, chain
import subprocess
from argparse import ArgumentParser

def get_expansion_info(file: str):
    pk = Pickaxe()
    pk.load_pickled_pickaxe(file)
    columns = [
        "starter_name",
        "starter_smiles",
        "starter_id",
        "target_name",
        "target_smiles",
        "target_id"
    ]

    right_order = [
        "starter_name",
        "target_name",
        "filename",
        "generations",
        "processed",
        "starter_smiles",
        "target_smiles",
        "starter_id",
        "target_id"
    ]

    generations = 0
    starters = set() # (name, smi, id)
    for v in pk.compounds.values():
        if v["Type"].startswith("Start"):
            starters.add((v["ID"], v["SMILES"], v["_id"]))

        if v["Generation"] > generations:
            generations = v["Generation"]
    
    targets = {(elt["ID"], elt["SMILES"], elt["_id"]) for elt in  pk.targets.values()}
    new_rows = [list(chain(*elt)) for elt in product(starters, targets,)]
    new_df = pd.DataFrame(data=new_rows, columns=columns)
    new_df["filename"] = file.stem
    new_df['generations'] = generations
    new_df["processed"] = False
    new_df = new_df[right_order]

    return new_df

def track_new_expansions(tracker: pd.DataFrame):
    tracked_filenames = set(tracker['filename'].to_list())

    for file in filepaths['raw_expansions'].iterdir():
        if file.stem not in tracked_filenames:
            new_rows = get_expansion_info(file)
            tracker = pd.concat((tracker, new_rows), axis=0)

    tracker.reset_index(drop=True, inplace=True)
    
    return tracker


def main(args):
    tracker = pd.read_csv(
        filepath_or_buffer=filepaths['artifacts'] / "expansion_tracking.csv",
        sep=',',
        index_col=0
    )
    tracker = track_new_expansions(tracker)
    tracker.to_csv(filepaths['artifacts'] / "expansion_tracking.csv", sep=',')

    gb = tracker.groupby("filename")
    for fn, group in gb:
        if not all(tracker.loc[group.index, "processed"]):
            gen = tracker.loc[group.index, "generations"].to_list()[0]

            try:
                if args.do_thermo:
                    subprocess.run(["python", "./scripts/process_expansion.py", fn, str(gen), "--do_thermo"], check=True)
                else:
                    subprocess.run(["python", "./scripts/process_expansion.py", fn, str(gen)], check=True)
            except subprocess.CalledProcessError as e:
                print(e)
            else:                
                tracker.loc[group.index, "processed"] = True
                tracker.to_csv(filepaths['artifacts'] / "expansion_tracking.csv", sep=',')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("--do_thermo", action="store_true", help="Does thermo calculations if provided")
    args = parser.parse_args()
    main(args)