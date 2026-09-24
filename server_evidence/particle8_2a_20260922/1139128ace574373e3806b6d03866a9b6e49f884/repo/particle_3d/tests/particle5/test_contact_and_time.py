import numpy as np
import pytest
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle3_cases import plane_triangle
from particle_3d.particle5_cases import contact_cases,synthetic_motion
from particle_3d.particle5_motion import simulate_resistance,resistance_trial
from particle_3d.particle4_motion import initial_world
from particle_3d.physical_time_refinement import refine_interval,PhysicalTimeRefinementError
from particle_3d.resistance_solver import ResistanceSystemIllConditioned


def test_resistance_metric_contact(mu):
    for row in contact_cases(mu):
        tol=row['p5_audit']['contact_kkt']['velocity_budget_m_s']
        assert min(row['p4_normal_speeds'])>=-tol and min(row['p5_normal_speeds'])>=-tol
        assert row['p5_resistance_objective']<row['p4_resistance_objective']
        assert row['velocity_difference']>1e-8
        assert row['p5_audit']['multiplier_role']=='KINEMATIC_CONSTRAINT_MULTIPLIER'


def test_resistance_contact_no_penetration_simultaneous_wall_pair(mu):
    shapes={1:Sphere([0,0,1e-6],1e-6),2:Sphere([0,0,3e-6],1e-6)}
    wall=WallGeometry([plane_triangle()]);provider=lambda i,s,t:(np.array([0,0,-i*1e-6]),np.zeros(3))
    end,rows,ledger=simulate_resistance(shapes,provider,mu,.1,.5,wall=wall)
    assert end.time_s==.5
    assert all(min(g['gap_m']+g['roundoff_budget_m'] for g in r['pair_gaps'])>=0 for r in rows)
    assert all(p['wall_gap_m']>=-p['wall_roundoff_m'] for r in rows for p in r['particles'])
    assert rows[-1]['projection']['contact_count']>=2


@pytest.mark.parametrize('kind',['wall','pair'])
def test_lubrication_physical_time_and_dt_comparison(kind,mu):
    final=[]
    for divisor in [1,2,4]:
        summary,states,ledger=synthetic_motion(mu,kind,divisor)
        assert summary['time_coverage_error_s']<1e-13
        assert summary['minimum_gap_m']>0
        assert summary['contact_events']==0
        for left,right in zip(ledger,ledger[1:]):assert abs(left['t1_s']-right['t0_s'])<1e-14
        assert all(r['projection']['dissipation']['total']>=0 for r in states[1:])
        final.append(summary['minimum_gap_m'])
    assert abs(final[2]-final[1])<abs(final[1]-final[0])


def test_unsafe_trial_splits_actual_physical_time(mu):
    shapes={1:Sphere([0,0,1.001e-6],1e-6)};wall=WallGeometry([plane_triangle()])
    provider=lambda i,s,t:(np.array([0,0,-1e-6]),np.zeros(3))
    end,rows,ledger=simulate_resistance(shapes,provider,mu,2.,2.,wall=wall)
    assert end.time_s==2. and any(l['depth']>0 for l in ledger)
    assert abs(sum(l['dt_s'] for l in ledger)-2.)<1e-14
    assert min(p['wall_gap_m']+p['wall_roundoff_m'] for r in rows for p in r['particles'])>=0


def test_resistance_solver_failure_does_not_consume_time(mu,monkeypatch):
    calls=[]
    def failed(*a,**k):
        calls.append(1);raise ResistanceSystemIllConditioned(dict(reason='TEST_NUMERICAL_FAILURE'))
    monkeypatch.setattr('particle_3d.particle5_motion.solve_resistance',failed)
    shape=Sphere([0,0,0],1e-6);provider=lambda i,s,t:(np.ones(3),np.zeros(3))
    old=initial_world({1:shape},provider);ledger=[];accepted=[]
    with pytest.raises(PhysicalTimeRefinementError) as exc:
        refine_interval(old,1.,resistance_trial(provider,mu,on_accept=accepted.append),max_depth=3,ledger=ledger)
    assert len(calls)==4 and not ledger and not accepted and old.time_s==0
    np.testing.assert_array_equal(old.shapes[1].center_m,[0,0,0])
    assert 'RESISTANCE_SYSTEM_ILL_CONDITIONED' in str(exc.value.__cause__)


def test_pair_tunneling_requests_real_subdivision(mu):
    a=1e-6;shapes={1:Sphere([-2*a,0,0],a),2:Sphere([2*a,0,0],a)}
    provider=lambda i,s,t:(np.array([4*a if i==1 else -4*a,0,0]),np.zeros(3))
    end,rows,ledger=simulate_resistance(shapes,provider,mu,1.,1.)
    assert end.time_s==1. and any(r['depth']>0 for r in ledger)
    assert end.shapes[1].center_m[0]<end.shapes[2].center_m[0]
    assert all(g['gap_m']>=-g['roundoff_budget_m'] for r in rows for g in r['pair_gaps'])
