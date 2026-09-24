from pathlib import Path
import csv
import json
import os
import sys
import pytest

PACKAGE = Path(__file__).resolve().parents[2]
REPO = PACKAGE.parent
sys.path.insert(0, str(PACKAGE / "src"))
from particle_3d.audit import read_frozen
from particle_3d.field import FrozenFEMField
from particle_3d.particle1_cases import uniform_case, rotation_case
from particle_3d.sonovue_adapter import sample_single_validation_size
from particle_3d.particle1_audit import VALIDATION_SEED


@pytest.fixture(scope="session")
def repo(): return REPO

@pytest.fixture(scope="session")
def sonovue_root(): return Path(os.environ.get("SONOVUE_ROOT", "/home/lzy/projects/sonovue_size_distribution_v0"))

@pytest.fixture(scope="session")
def size(sonovue_root): return sample_single_validation_size(sonovue_root, seed=VALIDATION_SEED)

@pytest.fixture(scope="session")
def uniform_rows(): return uniform_case()

@pytest.fixture(scope="session")
def rotation_rows(): return rotation_case()

@pytest.fixture(scope="session")
def p1_audited(): return read_frozen(REPO / "formal_3D_flow_solver/FEM_SimVascular")

@pytest.fixture(scope="session")
def p1_field(p1_audited): return FrozenFEMField.from_grids(p1_audited[1], p1_audited[2])

@pytest.fixture(scope="session")
def trajectory_cases():
    data = PACKAGE / "reports/particle1/data"
    result = []
    for index in range(3):
        with (data / f"05_real_trajectory_dt{index}.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        result.append((rows, json.loads((data / f"05_real_trajectory_dt{index}.json").read_text())))
    return result

@pytest.fixture(scope='session')
def trajectory_samples(trajectory_cases,p1_field):
    import numpy as np
    return [p1_field.sample_many(np.array([[float(r[f'{a}_m']) for a in 'xyz'] for r in rows])) for rows,_ in trajectory_cases]

@pytest.fixture(scope='session')
def classified_segments(trajectory_cases,p1_audited):
    import numpy as np
    from particle_3d.validation_boundary import ValidationBoundaryClassifier
    classifier=ValidationBoundaryClassifier(p1_audited[3])
    return [[classifier.first_event(np.array([float(a[f'{axis}_m']) for axis in 'xyz']),
                                     np.array([float(b[f'{axis}_m']) for axis in 'xyz']))
             for a,b in zip(rows[:-1],rows[1:])] for rows,_ in trajectory_cases]
