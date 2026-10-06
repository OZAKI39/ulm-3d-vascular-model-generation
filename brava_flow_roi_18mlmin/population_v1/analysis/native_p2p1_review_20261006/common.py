"""Isolated native higher-order review. Never writes production inputs."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]
PREVIOUS = HERE.parent / 'stabilization_local_conservation_review_20261006'
RETAINED = BASE / 'gpu_solver_fix/source'
SOURCE = RETAINED / 'Code/Source/solver'

def read(p):
    return json.loads(Path(p).read_text())

def put(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()
