from pathlib import Path
import os
import sys
import pytest

PACKAGE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PACKAGE / "src"))
from particle_3d.audit import read_frozen
from particle_3d.field import FrozenFEMField
from particle_3d.validation_cases import affine_case, real_node_case, shared_faces, shared_face_case, classification_case, real_flow_case


@pytest.fixture(scope="session")
def frozen_root():
    return Path(os.environ.get("PARTICLE0_FEM_ROOT", PACKAGE.parent / "formal_3D_flow_solver/FEM_SimVascular"))


@pytest.fixture(scope="session")
def audited(frozen_root):
    return read_frozen(frozen_root)


@pytest.fixture(scope="session")
def real_field(audited):
    return FrozenFEMField.from_grids(audited[1], audited[2])


@pytest.fixture(scope="session")
def affine():
    return affine_case()


@pytest.fixture(scope="session")
def node_rows(real_field, audited):
    return real_node_case(real_field, audited[3])


@pytest.fixture(scope="session")
def faces(real_field):
    return shared_faces(real_field)


@pytest.fixture(scope="session")
def face_case(real_field, faces):
    return shared_face_case(real_field, faces)


@pytest.fixture(scope="session")
def classifications(real_field, audited):
    return classification_case(real_field, audited[3])


@pytest.fixture(scope="session")
def random_real(real_field):
    return real_flow_case(real_field)
