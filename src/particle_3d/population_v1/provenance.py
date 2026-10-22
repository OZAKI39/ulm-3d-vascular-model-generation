"""Canonical provenance, explicit real-flow inputs, and reproducible records."""
from pathlib import Path
import hashlib
import json
import numpy as np


def canonical(value):
    def convert(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    return (json.dumps(value, sort_keys=True, separators=(',', ':'),
                       allow_nan=False, default=convert) + '\n').encode()


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def content_sha(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.pending')
    temp.write_bytes(canonical(value))
    temp.replace(path)
