from pathlib import Path
import csv
import json
import sys
import numpy as np
import pytest

PACKAGE=Path(__file__).resolve().parents[2]
REPO=PACKAGE.parent
DATA=PACKAGE/"reports/particle2/data"
sys.path.insert(0,str(PACKAGE/"src"))
from particle_3d.rbc_distribution import sample_rbc_geometries,load_contract,quantile_indices
from particle_3d.rbc import RBCGeometry


@pytest.fixture(scope="session")
def p2_repo():return REPO

@pytest.fixture(scope="session")
def p2_data():return DATA

@pytest.fixture(scope="session")
def contract():return load_contract()

@pytest.fixture(scope="session")
def population():return sample_rbc_geometries(100000,2026092002)

@pytest.fixture(scope="session")
def geometries(population):return [RBCGeometry.from_population(population,i) for i in quantile_indices(population.samples)]

@pytest.fixture(scope="session")
def static_case_result(geometries):
    from particle_3d.particle2_cases import static_case
    return static_case(geometries)

@pytest.fixture(scope="session")
def rotation_case_result(geometries):
    from particle_3d.particle2_cases import rigid_rotation_case
    return rigid_rotation_case(geometries)

@pytest.fixture(scope="session")
def p2_audited():
    from particle_3d.audit import read_frozen
    return read_frozen(REPO/"formal_3D_flow_solver/FEM_SimVascular")

@pytest.fixture(scope="session")
def p2_field(p2_audited):
    from particle_3d.field import FrozenFEMField
    return FrozenFEMField.from_grids(p2_audited[1],p2_audited[2])


def read_trajectory(path):
    with path.open(newline="") as stream:raw=list(csv.DictReader(stream))
    return [{k:v if k in ["timestep_role","boundary_event"] else v=="True" if k=="inside_lumen" else float(v)
             for k,v in row.items()} for row in raw]


@pytest.fixture(scope="session")
def real_cases():
    return [(read_trajectory(DATA/f"08_real_g{g}_dt{d}.csv"),json.loads((DATA/f"08_real_g{g}_dt{d}.json").read_text()))
            for g in range(5) for d in range(3)]

@pytest.fixture(scope="session")
def center_samples(real_cases,p2_field):
    # Shapes do not affect V=u; verify all 5 centers coincide in tests, sample
    # each unique center once using unchanged P0, without a new field cache.
    return [p2_field.sample_many(np.array([[row[f"{a}_m"] for a in "xyz"] for row in real_cases[d][0]])) for d in range(3)]
