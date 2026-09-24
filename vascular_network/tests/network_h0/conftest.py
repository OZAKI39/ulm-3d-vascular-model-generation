from pathlib import Path
import json
import pytest
from network_1d0d.idealized_h0 import load_h0, solve_operating_point

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'reports/a_network_1d0d_boundary_v2_idealized/data'

@pytest.fixture(scope='session')
def domain():
    return load_h0(ROOT/'reports/a_network_1d0d_boundary_v1')

@pytest.fixture(scope='session')
def baseline(domain):
    return solve_operating_point(domain)

@pytest.fixture(scope='session')
def bc():
    return json.loads((DATA/'roi_fixed_pressure_bc_H0.json').read_text())
