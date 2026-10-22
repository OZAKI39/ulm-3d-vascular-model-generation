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


def test_shared_clock_tiny_steps_only_stop_limiting_particles():
    torch=pytest.importorskip('torch')
    if not torch.cuda.is_available():pytest.skip('Requires CUDA')
    from particle_3d.population_v1.gpu_safety import GPUSafetyMonitor
    from particle_3d.population_v1.safety import CPUSafetyMonitor,PopulationSafetyEvent
    env=shear_environment();monitor=GPUSafetyMonitor(env,.0005,'cuda')
    x=torch.tensor([[0.,0.,2e-6],[.0001,0.,2e-6]],device='cuda',dtype=torch.float64)
    for row in x:monitor.insert(row,0.)
    for k in range(127):
        assert not bool(monitor.provider().any())
        v=torch.ones_like(x)*.01;end=x+1e-12*v
        monitor.accepted(x,end,v,(k+1)*1e-12,torch.ones(2,device='cuda',dtype=torch.bool),torch.tensor([True,False],device='cuda'))
        x=end
    reasons=monitor.provider().cpu().tolist()
    assert reasons[0] in (1,2) and reasons[1]==0,reasons
    cpu=CPUSafetyMonitor(env,.0005)
    from collections import deque
    for i in [1,2]:
        cpu.insert(particle(env,i,[0,0,2e-6]),0.)
        r=cpu.rows[i];r['calls']=127;r['count']=128
        r['history']=deque([(k*1e-12,np.array([k*1e-14,0,2e-6])) for k in range(64)],maxlen=64)
        r['limiter']=deque([i==1]*64,maxlen=64)
    with pytest.raises(PopulationSafetyEvent):cpu.provider(1)
    cpu.provider(2)
    save('shared_clock_guard_regression',dict(PASS=True,GPU_reasons=reasons,
        unchanged_thresholds=True,scope='Actual handoff owner only; unrelated synchronized particle survives'))


def test_rejected_interval_tags_only_actual_pair_owners():
    torch=pytest.importorskip('torch')
    if not torch.cuda.is_available():pytest.skip('Requires CUDA')
    from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
    from particle_3d.physical_time_refinement import TrialNeedsSubdivision
    env=shear_environment();engine=GPUPopulationEngine(env,pairs=True)
    for pid,x in [(1,[0,0,2e-6]),(2,[5e-6,0,2e-6]),(3,[.001,0,2e-6])]:
        p=particle(env,pid,x);engine.insert_accepted(p,event(p))
    velocity=torch.tensor([[.01,0,0],[-.01,0,0],[0,0,0]],device='cuda',dtype=torch.float64)
    with pytest.raises(TrialNeedsSubdivision):engine._pair_certificate(velocity,.0005,torch.tensor([[0,1]],device='cuda'))
    assert engine.refinement_owners=={1,2}
    save('rejected_trial_owner_regression',dict(PASS=True,responsible_particle_ids=sorted(engine.refinement_owners),unrelated_particle_id=3))
