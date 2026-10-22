"""Parity against the actual existing P9-A.1 stepper, not a substitute formula."""
from pathlib import Path
from types import SimpleNamespace
import os
import numpy as np
import pyvista as pv
from particle_3d.field import FrozenFEMField
from particle_3d.wall_geometry import WallGeometry
from particle_3d.validation_boundary import ValidationBoundaryClassifier
from particle_3d.particle_shapes import Sphere
from particle_3d.lammps_state import BridgeParticle
from particle_3d.lammps_neighbors import ValidationNeighborPolicy
from particle_3d.particle9a_motion import Particle9AStepper
from particle_3d.population_v1.engine import PopulationEngine
from particle_3d.population_v1.provenance import write_json


def shear_environment(outlet_x=None):
    points = np.array([[-2,-2,-1],[6,-2,-1],[-2,6,-1],[-2,-2,7]],float)*1e-3
    velocities = np.zeros((4,3)); velocities[:,0] = 1000*points[:,2]
    field = FrozenFEMField(points, np.array([[0,1,2,3]]), velocities, np.zeros(4))
    corners = np.array([[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0]],float)*1e-3
    wall = WallGeometry(corners[np.array([[0,2,1],[0,3,2]])],
                        provenance={'role':'SYNTHETIC_PLANAR_VALIDATION_ONLY'})
    classifier = None
    if outlet_x is not None:
        c = corners.copy(); c[:,2] = c[:,0]; c[:,0] = outlet_x
        s = pv.PolyData(c, np.array([[3,0,1,2],[3,0,2,3]]).ravel())
        classifier = ValidationBoundaryClassifier({'OUTLET_01':s})
    return SimpleNamespace(field=field, wall=wall, mu=.00345312, classifier=classifier)


def particle(env, pid, center, radius=1e-6):
    sample = env.field.sample(center)
    return BridgeParticle.from_shape(pid, Sphere(center,radius),
                                     velocity=sample.velocity_m_s,omega=.5*sample.vorticity_s_inv)


def event(p, birth=0.):
    return dict(particle_id=p.particle_id,source_event_id=p.particle_id,
                realization_id=0,birth_time_s=birth,clearance_m=0.)


def reference(env, p):
    def provider(i,s,t):
        sample=env.field.sample(s.center_m)
        return sample.velocity_m_s,.5*sample.vorticity_s_inv
    return Particle9AStepper([p],ValidationNeighborPolicy(20e-6,.5e-6,'PARITY_ONLY'),
        provider,env.mu,wall=env.wall,boundary_classifier=env.classifier,
        gradient_provider=lambda x:env.field.sample(x).velocity_gradient_s_inv)


def compare(env, particles, steps, dt=.0005):
    refs=[reference(env,p) for p in particles]
    engine=PopulationEngine(env,dt=dt,pairs=False)
    for p in particles:engine.insert_accepted(p,event(p))
    errors=dict(position_m=0.,velocity_m_s=0.,omega_s_inv=0.,wall_gap_m=0.)
    for step in range(1,steps+1):
        t=step*dt
        for ref in refs:
            if not ref.boundary_event.startswith('OUTLET_'):ref.step_to(t)
        engine.advance(t)
        for ref in refs:
            p=ref.read()[0]
            snapshots=[r for r in engine.snapshots if p.particle_id in r['particle_id']]
            row=snapshots[-1];k=int(np.flatnonzero(row['particle_id']==p.particle_id)[0])
            for key,a,b in [('position_m',row['x'][k],p.position),('velocity_m_s',row['v'][k],p.velocity),
                            ('omega_s_inv',row['omega'][k],p.omega),('wall_gap_m',row['wall_gap'][k],ref.samples[-1][14])]:
                errors[key]=max(errors[key],float(np.max(np.abs(np.asarray(a)-b))))
            if ref.boundary_event.startswith('OUTLET_'):
                terminal=next(r for r in engine.terminals if r['particle_id']==p.particle_id)
                assert terminal['outlet']==ref.boundary_event
                assert abs(terminal['residence_time_s']-ref.time_s)<1e-12
        if engine.pool.n==0:break
    assert errors['position_m']<max(env.wall.roundoff_m*8,1e-18),errors
    assert errors['velocity_m_s']<1e-11,errors
    assert errors['omega_s_inv']<1e-7,errors
    assert errors['wall_gap_m']<max(env.wall.roundoff_m*8,1e-18),errors
    return dict(PASS=True,errors=errors,steps=step,terminal_outlets=[r['outlet'] for r in engine.terminals],
                expected_outlets=[r.boundary_event for r in refs],active_count=engine.pool.n,
                flow_role=getattr(env,'flow_sha256','SYNTHETIC_PLANAR_VALIDATION_ONLY'))


def save(name, result):
    folder=Path(os.environ.get('POPULATION_VALIDATION_DIR','/tmp/population_v1_validation'))
    write_json(folder/(name+'.json'),result)


def test_single_particle_nearwall_parity():
    env=shear_environment(2e-6)
    p=particle(env,1,[0.,0.,1.01e-6])
    result=compare(env,[p],50)
    assert result['terminal_outlets']==['OUTLET_01']
    save('phase2_single_particle',result)


def test_multi_particle_noninteraction_parity():
    env=shear_environment()
    particles=[particle(env,1,[-1e-4,0,1.01e-6]),particle(env,2,[1e-4,0,1.005e-6])]
    result=compare(env,particles,20)
    save('phase2_noninteraction',result)


def test_so_a_and_exact_birth_time():
    env=shear_environment();engine=PopulationEngine(env)
    engine.advance(.000137)
    p=particle(env,8,[0,0,2e-6]);engine.insert_accepted(p,event(p,.000137))
    assert engine.pool.birth_time[0]==.000137
    engine.advance(.0005)
    assert engine.pool.age[0]==.0005-.000137
    assert engine.pool.x.dtype==np.float64
    assert engine.accounting()['conservation_pass']
    save('phase2_clock',dict(PASS=True,birth_time_s=.000137,accepted_time_s=engine.time))
