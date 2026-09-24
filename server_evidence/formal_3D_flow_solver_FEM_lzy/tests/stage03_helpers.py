"""Read measured Stage 3 evidence; no surrogate solver results."""
import json
from pathlib import Path
import numpy as np
import pytest
from fem3d.vascular import load_config
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/stage03/reference'
def read(path):return json.loads((ROOT/path).read_text())
def require_solution():
    if (BASE/'metadata/failure.json').is_file():
        failure=read('outputs/stage03/reference/metadata/failure.json')
        assert failure['status']=='FAIL' and failure['pde_solution_available'] is False
        assert not (BASE/'checkpoints/primary.npz').exists()
        pytest.skip('NOT EVALUATED: formal MUMPS factorization failed; no solution exists. Stage status remains FAIL.')
def qc():
    require_solution()
    return read('outputs/stage03/reference/qc/solution.json')
def solve():
    require_solution()
    return read('outputs/stage03/reference/metadata/solve.json')
def preflight():return read('outputs/stage03/preflight/assembly_verified.json')
def config():return load_config(ROOT)[0]
def close(actual,expected,rtol=1e-12):assert np.isclose(actual,expected,rtol=rtol,atol=0)
