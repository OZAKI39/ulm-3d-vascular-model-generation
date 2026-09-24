import numpy as np
import pytest
from particle_3d.particle9a_diagnostics import Observer, wrap_nominal_steps
from particle_3d.particle81_simulation import SavedTrajectoryStepper
from particle_3d.particle9a_motion import Particle9AStepper
from particle_3d.particle65_cases import synthetic, records_for, MU
from particle_3d.lammps_neighbors import ValidationNeighborPolicy

@pytest.mark.parametrize('model', ['P65', 'P9A'])
def test_logging_bitwise_parity_with_rejected_trials(model):
    shapes, _, wall, _ = synthetic('wall')
    def provider(i, s, t): return np.array([1e-5, 0., -1e-4]), np.array([0., 25., .3])
    gradient = np.zeros((3, 3)); gradient[0, 2] = 50.
    cls = SavedTrajectoryStepper if model == 'P65' else Particle9AStepper
    kw = {} if model == 'P65' else dict(gradient_provider=lambda x:gradient)
    def create(): return cls(records_for(shapes, provider), ValidationNeighborPolicy(20e-6, .5e-6, 'DIAGNOSTIC_TEST'), provider, MU, wall=wall, **kw)
    off, on = create(), create(); observer=Observer()
    assert wrap_nominal_steps(off, Observer()) is off
    wrap_nominal_steps(on, observer, enabled=True)
    for t in [.001, .004, .016, .032]:
        off.step_to(t); on.step_to(t)
        a, b = off.read()[0], on.read()[0]
        for name in ('position', 'velocity', 'omega', 'q'):
            assert np.array_equal(getattr(a, name), getattr(b, name))
    assert np.array_equal(np.array(off.samples), np.array(on.samples))
    assert on.rejected_count == off.rejected_count > 0
    assert any(not r['trial_accepted'] for r in observer.trials)
    assert max(r['subdivision_depth'] for r in observer.trials) > 0

@pytest.mark.parametrize('model', ['P65', 'P9A'])
def test_real_field_observation_bitwise_parity(model):
    """Exercise real FEM/geometry/closure observations, not only the hooks."""
    from pathlib import Path
    import json
    from copy import deepcopy
    from particle_3d.particle81_simulation import environment
    from particle_3d.injection_admission import FiniteSizeAdmission
    root=Path(__file__).resolve().parents[2]
    events=json.loads((root/'reports/particle9a_2mmps_diagnosis/data/original_smoke_audit.json').read_text())['events']
    events=list(events.values()) if isinstance(events,dict) else events
    event=next(e for e in events if e['particle_id']==12);env=environment()
    checker=FiniteSizeAdmission(None,None,wall=env.wall,field=env.field)
    def provider(i,s,t):
        f=env.field.sample(s.center_m);return f.velocity_m_s,.5*f.vorticity_s_inv
    cls=SavedTrajectoryStepper if model=='P65' else Particle9AStepper
    kw={} if model=='P65' else dict(gradient_provider=lambda x:env.field.sample(x).velocity_gradient_s_inv)
    def create():
        p=checker.check(deepcopy(event),np.array(event['birth_center_m']),{})[0]
        return cls([p],ValidationNeighborPolicy(20e-6,.5e-6,'DIAG_REAL_PARITY'),provider,env.mu,
                   wall=env.wall,boundary_classifier=env.classifier,**kw)
    off,on=create(),create();observer=Observer(env,event,model)
    wrap_nominal_steps(on,observer,enabled=True)
    for k in range(1,25):off.step_to(k*.00025);on.step_to(k*.00025)
    assert np.array_equal(off.samples,on.samples)
    assert len(observer.states)==len(on.samples)
    assert all('FEM_gradient_3x3' in r for r in observer.states)
