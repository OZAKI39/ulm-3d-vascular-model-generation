"""Captured SV1.3L evidence and synthetic rejection cases; no solver launches."""
import copy,json,math
from pathlib import Path
import pytest
from sv_validation.provenance import sha256
from sv_validation.sv13l import *
ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'reports/sv1_3l'
def load(name):return json.loads((R/(name+'.json')).read_text())
def actual(name):
    d=load(name)
    if d['status']=='NOT_RUN':
        assert d['executed'] is False and d['measurements'] is None
        pytest.skip(d['reason'])
    assert d['status']=='PASS',d.get('reason',d)
    return d
def policy():return json.loads((ROOT/'configs/sv1_3l/policy.json').read_text())
def benchmark_sample(t):return dict(wall_time_s=t,profiling=False,steps=20,initial_state='t=0',science_pass=True)
def proof_fixture():
    return dict(executed=True,steps=20,linear_failures=0,nonlinear_failures=0,
      velocity_finite=True,pressure_finite=True,Qin=3e-15,Qout={'OUTLET_01':1e-15,'OUTLET_02':1e-15,'OUTLET_03':1e-15},mass_error=0,wall_noslip=True,reload_pass=True)

