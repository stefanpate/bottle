from equilibrator_cache import ZenodoSettings
from hydra import initialize, compose
from pathlib import Path

with initialize(version_base=None, config_path="../conf/filepaths"):
    fps = compose(config_name="filepaths")

def get_cached_filepath(settings: ZenodoSettings) -> Path:
    cache_dir = Path(fps.equilibrator_cache)
    filepath = cache_dir / settings.filename

    if not filepath.exists():
        raise FileNotFoundError(f"Cached file not found: {filepath}. Ensure that equilibrator cache is saved to {cache_dir} and that conf/filepaths/filepaths.yaml is correct.")

    return cache_dir / settings.filename