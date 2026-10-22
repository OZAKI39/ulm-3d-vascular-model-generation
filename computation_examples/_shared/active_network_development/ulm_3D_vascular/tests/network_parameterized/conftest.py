from pathlib import Path
import pytest
from network_1d0d.parameterized_hydraulics import load_parameterized_domain, load_geometry_cache

ROOT=Path(__file__).resolve().parents[2]
GRAPH=ROOT/'reports/a_network_1d0d_boundary_v2_idealized/data/analysis_A_H0_graph_si.npz'
PORTS=ROOT/'reports/a_network_1d0d_boundary_v1/data/roi_ports_in_a.json'

@pytest.fixture(scope='session')
def geometry_domain():
    return load_parameterized_domain(GRAPH,PORTS)

@pytest.fixture(scope='session')
def geometry_cache():
    return load_geometry_cache(GRAPH,PORTS)
