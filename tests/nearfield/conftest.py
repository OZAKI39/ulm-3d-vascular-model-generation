import pytest
@pytest.fixture(scope='session')
def policy():
    from mouse_vessel_flow.particles.nearfield_regularization import NearFieldRegularizationV1
    return NearFieldRegularizationV1()
@pytest.fixture(scope='session')
def scans():
    from mouse_vessel_flow.particles.particle65_validation import static_scans
    return static_scans()
@pytest.fixture(scope='session')
def wall_run():
    from mouse_vessel_flow.particles.particle65_cases import run_case
    return run_case('wall',.004,.032)
@pytest.fixture(scope='session')
def pair_run():
    from mouse_vessel_flow.particles.particle65_cases import run_case
    return run_case('pair',.004,.032)
@pytest.fixture(scope='session')
def bridge_result(tmp_path_factory):
    from mouse_vessel_flow.particles.particle65_validation import bridge_case
    return bridge_case(tmp_path_factory.mktemp('v1_binary')/'checkpoint')
