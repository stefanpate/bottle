import json
import orjson
from pathlib import Path
from typing import TypeVar

T = TypeVar('T')

def save_json(data, save_to):
    with open(save_to, 'w') as f:
        json.dump(data, f)

def load_json(path):
    json_str = Path(path).read_text('utf-8')
    data = orjson.loads(json_str)
    return data