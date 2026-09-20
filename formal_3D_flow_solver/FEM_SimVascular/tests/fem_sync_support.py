import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT/name).read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
