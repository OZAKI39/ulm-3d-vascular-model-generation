from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from particle_3d.injection_method_c import TruncatedSonoVue,MethodCSource
from particle_3d.inlet_flux import InletFluxSampler
from particle_3d.particle7_cases import channel,SONOVUE

@pytest.fixture(scope='session')
def distribution():return TruncatedSonoVue(SONOVUE)

from method_c_test_helpers import source

@pytest.fixture
def synthetic(distribution):return source(distribution)
