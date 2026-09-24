from pathlib import Path
import pytest
from particle_3d.routing_stationary_audit import read
ROOT=Path(__file__).resolve().parents[3]
REPORT=ROOT/'particle_3d/reports/particle9a1_routing_stationary_audit'
@pytest.fixture(scope='session')
def audit():return REPORT
@pytest.fixture(scope='session')
def root():return ROOT
@pytest.fixture(scope='session')
def events():return read(ROOT/'particle_3d/outputs/particle9a1_2mmps/admission/birth_ledger.json')['events']
@pytest.fixture(scope='session')
def candidates():return read(ROOT/'particle_3d/outputs/particle9a_2mmps/admission/all_scheduled_admission_records.json')
