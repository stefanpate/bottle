import argparse
from pathlib import Path
import yaml

def main(args):
    project_dir = Path(__file__).parent.parent
    in_project_paths = ["scripts", "logs", "artifacts"]
    
    filepaths = {
        'dirs': {
            name: str(project_dir / name) for name in in_project_paths
        },
        'subdirs': {
            'results': [
                "processed_expansions",
                "raw_expansions",
            ],
            'artifacts': [
                "cofactors",
                "rules",
                "starters_targets",
                "imgs",
                "operator_mapping",
            ],
            'data': [
                "sprhea",
            ],
        },
    }

    # Add user selected top-level directories
    filepaths['dirs']['data'] = args.data
    filepaths['dirs']['results'] = args.results

    # Make sub directories for directories external to project directory
    for sub in filepaths['subdirs']['data']:
        if not (Path(args.data) / sub).exists():
            (Path(args.data) / sub).mkdir(parents=True, exist_ok=False)

    for sub in filepaths['subdirs']['results']:
        if not (Path(args.results) / sub).exists():
            (Path(args.results) / sub).mkdir(parents=True, exist_ok=False)

    # Make gitignored project subdirs
    if not Path(filepaths['dirs']['logs']).exists():
        Path(filepaths['dirs']['logs']).mkdir()

    # Write the filepaths to config.yaml
    config_path = Path(__file__).parent.parent / "config.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(filepaths, f, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Configs initializer.")
    parser.add_argument('data', type=str, help='Filepath to data directory')
    parser.add_argument('results', type=str, help='Filepath to results directory')

    args = parser.parse_args()
    
    main(args)