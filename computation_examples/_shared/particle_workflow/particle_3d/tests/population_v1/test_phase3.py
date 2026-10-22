import numpy as np
from types import SimpleNamespace
from particle_3d.field import FrozenFEMField
from particle_3d.particle_shapes import Sphere
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.resistance_solver import solve_resistance
from particle_3d.nearfield_handoff import handoff_constraints
from particle_3d.particle6_stepper import bind_query_dependency
from particle_3d.population_v1.pairs import candidates,sphere_gap,components,contract
from particle_3d.population_v1.components import sparse_assembly,solve_large
from particle_3d.population_v1.engine import PopulationEngine
from test_phase2 import particle,event,save,shear_environment,reference


def compression_environment():
    points=np.array([[-1,-1,-1],[3,-1,-1],[-1,3,-1],[-1,-1,3]],float)*1e-3
    # Incompressible affine convergence: dux/dx=-1000, duy/dy=1000.
    velocity=points*np.array([-1000.,1000.,0.])
    return SimpleNamespace(field=FrozenFEMField(points,np.array([[0,1,2,3]]),velocity,np.zeros(4)),
        wall=None,mu=.00345312,classifier=None)


def collide(dt,oblique=False):
    env=compression_environment();e=PopulationEngine(env,dt=dt,pairs=True)
    y=.15e-6 if oblique else 0.
    for pid,p in [(19,[-1.02e-6,y,0]),(81,[1.02e-6,-y,0])]:
        b=particle(env,pid,p);e.insert_accepted(b,event(b))
    minimum=float('inf');contacts=0;accepted=0
    class Recorder:
        def snapshot(self,row):pass
        def audit(self,world):
            nonlocal minimum,contacts,accepted
            minimum=min(minimum,sphere_gap(world.shapes[19],world.shapes[81],19,81).gap_m)
            contacts=max(contacts,world.projection['contact_count']);accepted+=1
    e.recorder=Recorder()
    e.step_to(.008)
    assert list(e.pool.particle_id[:2])==[19,81]
    assert minimum>=2e-9-2e-18
    assert e.accounting()['conservation_pass']
    return e.pool.x[:2].copy(),dict(dt_s=dt,minimum_pair_gap_m=minimum,
        maximum_contact_count=contacts,accepted_intervals=accepted,position_m=e.pool.x[:2].tolist())


def test_pair_collision_and_refinement():
    runs=[]
    for oblique in [False,True]:
        positions=[];rows=[]
        for dt in [.0005,.00025,.000125]:
            x,r=collide(dt,oblique);positions.append(x);rows.append(r)
        errors=[float(np.max(np.abs(positions[i]-positions[-1]))) for i in [0,1]]
        assert errors[1]<=errors[0]+1e-18
        if not oblique:assert rows[0]['maximum_contact_count']>0
        runs.append(dict(oblique=oblique,runs=rows,position_difference_to_finest_m=errors))
    save('phase3_collision',dict(PASS=True,fixtures=runs,contract=contract()))


def test_swept_broadphase_detects_fast_crossing():
    shapes={3:Sphere([-1e-4,0,0],1e-6),9:Sphere([1e-4,0,0],1e-6)}
    assert candidates(shapes)==[]
    assert candidates(shapes,velocities={3:np.array([.4,0,0]),9:np.array([-.4,0,0])},dt=.0005)==[(3,9)]


def test_large_sparse_chain_matches_original_dense_reference():
    policy=NearFieldRegularizationV1();n=40
    shapes={i:Sphere([i*(2e-6+2e-9),0,0],1e-6) for i in range(n)}
    free={i:np.array([.001*(n/2-i)/n,2e-4,0,0,0,0]) for i in shapes}
    system=sparse_assembly(shapes,free,.00345312,None,policy=policy)
    constraints,_=bind_query_dependency(handoff_constraints,dict(all_pairs=candidates,pair_gap=sphere_gap))(shapes,None,policy)
    result=solve_large(system,particles=shapes,constraints=constraints)
    ref=solve_resistance(system,particles=shapes,constraints=constraints)
    error=float(np.max(np.abs(result.velocity-ref.velocity)))
    assert error<1e-10,error
    assert result.record['max_constraint_violation_m_s']<=result.record['contact_kkt']['velocity_budget_m_s']
    save('phase3_sparse_chain',dict(PASS=True,particle_count=n,contacts=len(constraints),
        maximum_velocity_error_m_s=error,sparse_audit=result.record))


def test_components_preserve_isolated_wall_physics():
    env=shear_environment();engine=PopulationEngine(env,pairs=True)
    bubbles=[particle(env,3,[-1e-4,0,1.01e-6]),particle(env,90,[1e-4,0,1.005e-6])]
    refs=[reference(env,b) for b in bubbles]
    for b in bubbles:engine.insert_accepted(b,event(b))
    error=0.
    for k in range(1,11):
        t=k*.0005;engine.advance(t)
        for i,r in enumerate(refs):
            r.step_to(t);error=max(error,float(np.max(abs(engine.pool.x[i]-r.read()[0].position))))
    assert error<1e-18
    save('phase3_components',dict(PASS=True,isolated_position_error_m=error,
        components=components([3,90],[]),global_dense_matrix_created=False))
