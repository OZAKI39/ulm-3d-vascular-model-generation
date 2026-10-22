import os
from pathlib import Path
import numpy as np
import pytest
torch=pytest.importorskip('torch')
pytestmark=pytest.mark.skipif(not torch.cuda.is_available(),reason='Requires actual CUDA')
from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
from particle_3d.population_v1.engine import PopulationEngine
from particle_3d.population_v1.source import PerfusionSource,perfuse_until
from particle_3d.population_v1.provenance import write_json,canonical
from test_phase2 import shear_environment,particle,event
from test_phase3 import compression_environment
from test_phase4 import source_env


def save(name,row):write_json(Path(os.environ.get('POPULATION_VALIDATION_DIR','/tmp/population_v1_validation'))/(name+'.json'),row)


@pytest.mark.parametrize('fixture',['nearwall_outlet','separated','collision','oblique'])
def test_complete_cuda_integrator_parity(fixture):
    env=shear_environment(2e-6) if fixture=='nearwall_outlet' else compression_environment() if fixture in ('collision','oblique') else shear_environment()
    points=[[0.,0.,1.01e-6]] if fixture=='nearwall_outlet' else [[-1e-4,0,1.01e-6],[1e-4,0,1.005e-6]]
    if fixture in ('collision','oblique'):
        y=.15e-6 if fixture=='oblique' else 0.;points=[[-1.02e-6,y,0],[1.02e-6,-y,0]]
    cpu=PopulationEngine(env,pairs=True);gpu=GPUPopulationEngine(env,pairs=True)
    for k,x in enumerate(points):
        p=particle(env,k+1,x);cpu.insert_accepted(p,event(p));gpu.insert_accepted(p,event(p))
    errors=dict(position_m=0.,velocity_m_s=0.,omega_s_inv=0.,wall_gap_m=0.)
    for step in range(1,21):
        target=.0005*step;cpu.advance(target);gpu.advance(target)
        assert cpu.pool.n==gpu.pool.n
        if not cpu.pool.n:break
        for name,field in [('position_m','x'),('velocity_m_s','v'),('omega_s_inv','omega'),('wall_gap_m','wall_gap')]:
            a=getattr(cpu.pool,field)[:cpu.pool.n];b=getattr(gpu.pool,field)[:gpu.pool.n].cpu().numpy()
            errors[name]=max(errors[name],float(np.max(abs(a-b))))
    assert errors['position_m']<1e-12 and errors['velocity_m_s']<1e-10 and errors['omega_s_inv']<1e-4,errors
    assert [r['outlet'] for r in cpu.terminals]==[r['outlet'] for r in gpu.terminals]
    for a,b in zip(cpu.terminals,gpu.terminals):assert abs(a['exit_time_s']-b['exit_time_s'])<1e-10
    save('phase7_engine_'+fixture,dict(PASS=True,fixture=fixture,errors=errors,
        cpu_terminals=cpu.terminals,gpu_terminals=gpu.terminals,backend_counts=dict(gpu.backend_counts)))


def test_cuda_source_ledger_and_population_conservation():
    env=source_env();cpu=PopulationEngine(env,pairs=True);gpu=GPUPopulationEngine(env,pairs=True)
    a=PerfusionSource(env,concentration_m3=2e14,seed=4);b=PerfusionSource(env,concentration_m3=2e14,seed=4)
    perfuse_until(cpu,a,.003);perfuse_until(gpu,b,.003)
    assert canonical(cpu.source_ledger)==canonical(gpu.source_ledger)
    error=float(np.max(abs(cpu.pool.x[:cpu.pool.n]-gpu.pool.x[:gpu.pool.n].cpu().numpy())))
    assert error<1e-12
    save('phase7_source_gpu',dict(PASS=True,source_ledger_bit_identical=True,
        maximum_position_error_m=error,accounting=gpu.accounting(),backend_counts=dict(gpu.backend_counts)))
