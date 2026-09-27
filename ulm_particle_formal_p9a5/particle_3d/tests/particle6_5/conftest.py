from pathlib import Path
import sys,json
import pytest
PACKAGE=Path(__file__).resolve().parents[2];sys.path.insert(0,str(PACKAGE/'src'))
@pytest.fixture(scope='session')
def repo():return PACKAGE.parent
@pytest.fixture(scope='session')
def report():return PACKAGE/'reports/particle6_5'
@pytest.fixture(scope='session')
def policy():
    from particle_3d.nearfield_regularization import NearFieldRegularizationV1
    return NearFieldRegularizationV1()
@pytest.fixture(scope='session')
def scans():
    from particle_3d.particle65_validation import static_scans
    return static_scans()
@pytest.fixture(scope='session')
def wall_run():
    from particle_3d.particle65_cases import run_case
    return run_case('wall',.004,.032)
@pytest.fixture(scope='session')
def pair_run():
    from particle_3d.particle65_cases import run_case
    return run_case('pair',.004,.032)
@pytest.fixture(scope='session')
def bridge_result(tmp_path_factory):
    from particle_3d.particle65_validation import bridge_case
    return bridge_case(tmp_path_factory.mktemp('v1_binary')/'checkpoint')
@pytest.fixture(scope='session')
def real_saved(report):return json.loads((report/'data/07_real_v1_dt1.json').read_text())
