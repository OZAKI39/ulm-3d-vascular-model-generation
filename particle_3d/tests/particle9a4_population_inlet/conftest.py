from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from particle_3d.continuous_infusion import ContinuousInfusionSource
from particle_3d.injection_method_c import TruncatedSonoVue
from particle_3d.inlet_flux import InletFluxSampler
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.particle7_cases import channel
from particle_3d.population_inlet_p9a4 import load_new_environment,make_contract,make_source
ROOT=Path(__file__).resolve().parents[3]

class FixedSize:
    def __init__(self,d=1e-6):self.d=d;self.calls=0
    def sample(self,rng):self.calls+=1;return self.d,dict(uniform=float(rng.random()))

class CountSampler:
    Q_m3_s=1.
    def __init__(self):self.calls=0
    def sample(self,rng):self.calls+=1;return np.array([[rng.random(),0.,0.]]),np.array([0])

class RejectChecker:
    wall=None
    def __init__(self):self.calls=0
    def check(self,event,position,active):
        assert active=={};self.calls+=1
        return None,'WALL_REJECTED',{}

class HalfChecker(RejectChecker):
    def check(self,event,position,active):
        assert active=={};self.calls+=1
        return (SimpleNamespace(), 'ACCEPTED', {}) if position[0]<.5 else (None,'WALL_REJECTED',{})

@pytest.fixture
def simple_source():
    return ContinuousInfusionSource(sampler=CountSampler(),distribution=FixedSize(),checker=HalfChecker(),
        concentration_m3=3.,master_seed=1594,flow_sha256='SYNTHETIC',source_contract_sha256='fixed',synthetic=True)

@pytest.fixture(scope='session')
def distribution():return TruncatedSonoVue(ROOT/'sonovue_size_distribution_v0')

@pytest.fixture(scope='session')
def real_env():return load_new_environment(ROOT)

@pytest.fixture(scope='session')
def real_source(real_env):return make_source(real_env,make_contract(ROOT,real_env))

class PoiseuilleSampler:
    def __init__(self,R=2e-6):self.R=R;self.Q_m3_s=1.
    def sample(self,rng,n=1):
        # Full Poiseuille radial CDF G(r)=2(r/R)^2-(r/R)^4.
        r=self.R*np.sqrt(1-np.sqrt(1-rng.random(n)));theta=rng.uniform(0,2*np.pi,n)
        return np.column_stack((r*np.cos(theta),r*np.sin(theta),np.zeros(n))),np.zeros(n,dtype=int)

class TwoSizes:
    def sample(self,rng):
        u=float(rng.random());return (.6e-6 if u<.4 else 2.4e-6),dict(uniform=u)

class CircleCheck:
    wall=None
    R=2e-6
    h=2e-9
    def check(self,event,position,active):
        assert not active
        accepted=np.linalg.norm(position[:2])+event['radius_m']+self.h<=self.R
        return (SimpleNamespace() if accepted else None),'ACCEPTED' if accepted else 'WALL_REJECTED',{}

@pytest.fixture(scope='session')
def poiseuille():
    from particle_3d.population_inlet_p9a4 import generate
    source=ContinuousInfusionSource(sampler=PoiseuilleSampler(),distribution=TwoSizes(),checker=CircleCheck(),
        concentration_m3=5.,master_seed=2594,flow_sha256='ANALYTIC_POISEUILLE',source_contract_sha256='two-size',synthetic=True)
    return generate(source,40000,workers=1)
