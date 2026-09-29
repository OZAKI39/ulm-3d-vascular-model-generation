from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from particle_3d.field import FrozenFEMField
from particle_3d.wall_geometry import WallGeometry
from particle_3d.population_v1.engine import PopulationEngine
from particle_3d.population_v1.output import PopulationOutput,iter_states
from particle_3d.population_v1.provenance import canonical
from test_phase2 import particle,event,save,shear_environment


def corner_environment():
    points=np.array([[-1,-1,-1],[3,-1,-1],[-1,3,-1],[-1,-1,3]],float)*1e-3
    field=FrozenFEMField(points,np.array([[0,1,2,3]]),np.full((4,3),-.001),np.zeros(4))
    triangles=[]
    for axis in range(3):
        t=np.zeros((3,3));ij=[k for k in range(3) if k!=axis]
        t[:,ij]=np.array([[-1,-1],[1,-1],[0,1]])*1e-3
        if np.cross(t[1]-t[0],t[2]-t[0])[axis]>0:t=t[[0,2,1]]
        triangles.append(t)
    return SimpleNamespace(field=field,wall=WallGeometry(np.array(triangles)),mu=.00345312,classifier=None)


def run_corner(backend):
    env=corner_environment();engine=backend(env,pairs=True)
    p=particle(env,1,np.full(3,1e-6+2e-9));engine.insert_accepted(p,event(p))
    engine.step_to(.1)
    assert engine.pool.n==0
    terminal=engine.terminals[0]
    assert terminal['terminal_reason']=='MODEL_SUPPORTED_STATIONARY_TERMINATION'
    assert terminal['outlet'] is None and terminal['support_verified']
    assert len(terminal['support_evidence'])==32
    assert engine.counts['solver_failed_count']==0 and engine.accounting()['conservation_pass']
    return terminal


def test_original_stationary_support_terminates_cpu_population():
    t=run_corner(PopulationEngine)
    save('phase4_stationary_cpu',dict(PASS=True,terminal=t,active_pool_removed=True))


def test_original_stationary_support_terminates_gpu_population():
    torch=pytest.importorskip('torch')
    if not torch.cuda.is_available():pytest.skip('Requires CUDA')
    from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
    t=run_corner(GPUPopulationEngine)
    save('phase7_stationary_gpu',dict(PASS=True,terminal=t,active_pool_removed=True))


def test_full_accepted_time_output_roundtrip(tmp_path):
    env=shear_environment();writer=PopulationOutput(tmp_path/'run',dict(flow_sha256='EXPLICIT_SYNTHETIC'),chunk_rows=3)
    engine=PopulationEngine(env,pairs=True,recorder=writer)
    p=particle(env,9,[0,0,2e-6]);engine.insert_accepted(p,event(p))
    engine.step_to(.002);expected=engine.pool.x[:1].copy();manifest=writer.close(engine.accounting())
    states=list(iter_states(tmp_path/'run'))
    assert len(states)==5 and np.array_equal(states[-1]['x'],expected)
    assert states[-1]['physical_time_s']==.002 and states[-1]['age'][0]==.002
    assert len(manifest['chunks'])==2
