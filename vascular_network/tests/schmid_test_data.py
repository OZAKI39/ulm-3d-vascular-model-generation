from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np


def write_synthetic_schmid(directory: Path, *, equal_pressure_edge: bool = False) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    coordinates = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, -1.0, 0.0],
            [3.0, 0.0, 0.0],
        ]
    )
    pressure = [10.0, 8.0, 6.0, 6.0, 4.0]
    if equal_pressure_edge:
        pressure[1] = pressure[0]
    tuples = [(1, 0), (1, 2), (1, 3), (2, 4), (3, 4)]
    flows = [2.0, 1.0, 1.0, 1.0, 1.0]
    points = [np.vstack((coordinates[u], coordinates[v])) for u, v in tuples]
    lengths = [float(np.linalg.norm(coordinates[u] - coordinates[v])) for u, v in tuples]
    vertices = {
        "pressure": pressure,
        "coords": list(coordinates),
        "pBC": [10.0, None, None, None, 4.0],
    }
    edges = {
        "diameter": [2.0] * len(tuples),
        "tuple": tuples,
        "flow": flows,
        "httBC": [None] * len(tuples),
        "nkind": [4] * len(tuples),
        "length": lengths,
        "htt": [0.4] * len(tuples),
        "nRBC": [1.0] * len(tuples),
        "diameters": [np.asarray([2.0, 2.0])] * len(tuples),
        "points": points,
    }
    for name, payload in (("verticesDict.pkl", vertices), ("edgesDict.pkl", edges)):
        with (directory / name).open("wb") as stream:
            pickle.dump(payload, stream, protocol=2)
    return directory
