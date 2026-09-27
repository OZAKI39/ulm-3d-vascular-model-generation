from copy import deepcopy
from types import SimpleNamespace
import numpy as np
import pytest
from particle_3d.inlet_flux import InletFluxSampler
from particle_3d.particle7_cases import channel, SONOVUE
from particle_3d.sonovue_adapter import read_sonovue
from particle_3d.particle82a_admission import common_event
from particle_3d.particle82a_geometry import segment_distance, perimeter_edges


@pytest.fixture(scope='session')
def simple_context():
    sampler,wall,_,_=channel(width=3e-6,length=20e-6)
    distribution=read_sonovue(SONOVUE)[1]
    class Geometry:
        search_horizon=6e-6
        def distances(self,p):
            return wall.nearest_center_triangle(p)[1],float(1.5e-6-max(abs(p[0]),abs(p[1]))),0.
        def aperture_radius(self,p):return min(self.distances(p)[:2])
        def full_margin(self,p,r):return float(min(p[2]-r,wall.nearest_center_triangle(p)[1]-r-2e-9))
    field=SimpleNamespace(locate=lambda p:(0,None))
    return SimpleNamespace(env=SimpleNamespace(sampler=sampler,wall=wall,field=field),geometry=Geometry(),distribution=distribution)


@pytest.fixture
def common_inputs(simple_context):
    return [common_event(i,ctx=simple_context) for i in range(1,9)]
