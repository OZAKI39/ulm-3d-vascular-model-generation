"""Explicit physical settings and project-relative paths, loaded once per process."""

from pathlib import Path
import os
import tomllib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = Path(os.environ.get("CLOT_CONFIG", ROOT / "configs/model.toml")).resolve()
CONFIG = tomllib.loads(CONFIG_PATH.read_text())
INPUTS = ROOT / "inputs"
CACHE = ROOT / "results/cache"
NATIVE = ROOT / "native"
LIB = Path(
    os.environ.get("CLOT_LAMMPS_LIBRARY", NATIVE / "build/liblammps.so")
).resolve()
MATERIAL = CONFIG["material"]
DX = MATERIAL["nominal_spacing"]
V = DX**3
RHO = MATERIAL["density"]
K = MATERIAL["bulk_modulus"]
G = MATERIAL["shear_modulus"]
HORIZON = MATERIAL["horizon"]
S0 = np.array(MATERIAL["onset_by_type"])
SC = np.array(MATERIAL["complete_by_type"])
BASE = CONFIG["time"]["base_step"]
MIN_H = BASE / CONFIG["time"]["minimum_divisor"]
END, RISE, HOLD, OFF, COHORT = [
    CONFIG["time"][key] for key in ("end", "rise", "hold", "off", "cohort")
]
MU = CONFIG["loading"]["viscosity"]
H_EFF = CONFIG["loading"]["effective_thickness"]
BETA = MU / H_EFF
NOMINAL_SHEAR = CONFIG["loading"]["nominal_shear"]
LABEL = "MECHANICAL_SHEAR_SURROGATE"
HET_LABEL = "EXPLORATORY_UNCALIBRATED"
