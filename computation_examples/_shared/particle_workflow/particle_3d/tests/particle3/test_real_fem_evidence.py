import json
from pathlib import Path
import numpy as np
import pytest
from particle_3d.particle_shapes import Sphere,Capsule,Ellipsoid,EPS
from particle_3d.rbc_orientation import rotation_matrix
from particle_3d.wall_gap import wall_gap


def finite_tree(value):
    if isinstance(value,dict):return all(finite_tree(v) for v in value.values())
    if isinstance(value,list):return all(finite_tree(v) for v in value)
    return bool(np.isfinite(value)) if isinstance(value,(int,float)) else True


@pytest.mark.parametrize('dt_index',range(3))
def test_real_mb_wall_state_finite_and_unchanged_initialization(p3_repo,real_wall,dt_index):
    data=p3_repo/'particle_3d/reports/particle3/data';prefix=f'10_real_mb_dt{dt_index}'
    rows=json.loads((data/(prefix+'_states.json')).read_text());summary=json.loads((data/(prefix+'_summary.json')).read_text())
    initial=json.loads((p3_repo/'particle_3d/reports/particle1/data/05_real_initialization.json').read_text())
    size=json.loads((p3_repo/'particle_3d/reports/particle1/data/01_single_mb_size_provenance.json').read_text())
    assert rows and summary['status']=='OUTLET_02' and finite_tree(rows)
    assert rows[0]['center_m']==initial['initial_position_m']
    assert summary['validation_dt_s']==initial['validation_timesteps_s'][-1]/2**dt_index
    assert all(r['radius_m']==size['radius_m'] and r['shape_mode']=='SPHERE_MB' for r in rows)
    for r in rows:assert r['wall_gap_m']>=-r['roundoff_m']
    for i in [0,len(rows)//2,len(rows)-1]:
        r=rows[i];gap=wall_gap(Sphere(r['center_m'],r['radius_m']),real_wall)
        assert abs(gap.gap_m-r['wall_gap_m'])<=gap.roundoff_m
    for before,after in zip(rows,rows[1:]):
        exact=np.array(before['center_m'])+(after['time_s']-before['time_s'])*np.array(after['corrected_velocity_m_s'])
        assert np.linalg.norm(np.array(after['center_m'])-exact)<=8*EPS*np.linalg.norm(exact)


@pytest.mark.parametrize('geometry_index',range(5))
def test_real_rbc_wall_state_shape_mode_volume_area_and_axis(p3_repo,real_wall,geometry_index):
    data=p3_repo/'particle_3d/reports/particle3/data'
    for di in range(3):
        prefix=f'11_real_rbc_g{geometry_index}_dt{di}';rows=json.loads((data/(prefix+'_states.json')).read_text());summary=json.loads((data/(prefix+'_summary.json')).read_text())
        assert finite_tree(rows) and finite_tree(summary)
        assert summary['status'] in ['OUTLET_01','OUTLET_02','OUTLET_03','DEFORMATION_SURROGATE_INFEASIBLE','VALIDATION_HORIZON_REACHED']
        if not rows:
            assert summary['status']=='DEFORMATION_SURROGATE_INFEASIBLE'
            assert summary['minimum_gap_m'] is None and summary['last_accepted_time_s']==0
            assert summary['failure']['accepted_state'] is False
            continue
        for row in rows:
            assert row['wall_gap_m']>=-row['roundoff_m']
            assert row['shape_mode'] in ['FREE_OBLATE','CAPILLARY_DEFORMED']
            if row['shape_mode']=='CAPILLARY_DEFORMED':
                assert abs(row['capsule_volume_m3']/row['original_volume_m3']-1)<16*EPS
                assert row['area_ratio']<=1+512*EPS
                assert row['jeffery_orientation_interpreted'] is False
                velocity=np.array(row['free_velocity_m_s']);np.testing.assert_allclose(row['capsule_axis_world'],velocity/np.linalg.norm(velocity),rtol=0,atol=8*EPS)
        for before,after in zip(rows,rows[1:]):
            if before['shape_mode']=='CAPILLARY_DEFORMED':assert after['q']==before['q']
        for i in sorted(set([0,len(rows)//2,len(rows)-1])):
            row=rows[i]
            shape=Capsule(row['center_m'],row['capsule_axis_world'],row['R_cap_m'],row['L_cap_m']) if row['shape_mode']=='CAPILLARY_DEFORMED' else Ellipsoid(row['center_m'],[row[k] for k in ['a_m','b_m','c_m']],rotation_matrix(row['q']))
            gap=wall_gap(shape,real_wall);assert abs(gap.gap_m-row['wall_gap_m'])<=gap.roundoff_m


def test_all_real_intervals_cover_only_accepted_time(p3_repo):
    data=p3_repo/'particle_3d/reports/particle3/data';summaries=list(data.glob('*_real_*_summary.json'));assert len(summaries)==18
    for path in summaries:
        summary=json.loads(path.read_text());intervals=json.loads(path.with_name(path.name.replace('_summary','_intervals')).read_text())
        assert all(r['dt_s']>0 and r['accepted'] and r['proof'] for r in intervals)
        assert all(abs(a['t1_s']-b['t0_s'])<=8*EPS for a,b in zip(intervals,intervals[1:]))
        assert abs(sum(r['dt_s'] for r in intervals)-summary['last_accepted_time_s'])<=16*EPS
        assert not summary['failed_interval_consumed_time']


def test_upstream_orientation_does_not_claim_deformed_jeffery_convergence(p3_repo):
    comparison=json.loads((p3_repo/'particle_3d/reports/particle3/data/12_timestep_comparison.json').read_text())
    assert len(comparison['cases'])==18 and len(comparison['upstream_orientation'])==10
    assert comparison['production_particle_timestep_frozen'] is False
    for row in comparison['upstream_orientation']:
        if not row['upstream_matched_states']:
            assert row['max_upstream_axis_difference_deg'] is None
            assert row['status']=='NO_UPSTREAM_FREE_OBLATE_SEGMENT'
        assert row['particle2_orientation_convergence']=='NOT_ESTABLISHED'


def test_complete_cached_replay_matches_preserved_uncached_reference(p3_repo):
    report=p3_repo/'particle_3d/reports/particle3'
    reference=json.loads((report/'logs/uncached_real_rbc_g0_dt0_states.json').read_text())
    current=json.loads((report/'data/11_real_rbc_g0_dt0_states.json').read_text())
    assert len(reference)==len(current)
    for before,after in zip(reference,current):
        assert before['time_s']==after['time_s'] and before['shape_mode']==after['shape_mode']
        assert before['q']==after['q']
        assert np.linalg.norm(np.array(before['center_m'])-after['center_m'])<=after['roundoff_m']
        # Centers can differ at summation roundoff; compare geometry at the
        # original coordinate budget, not a newly enlarged wall tolerance.
        assert abs(before['R_cap_m']-after['R_cap_m'])<=after['roundoff_m']
        assert abs(before['wall_gap_m']-after['wall_gap_m'])<=after['roundoff_m']


def test_all_actual_contact_velocities_satisfy_predeclared_roundoff_bound(p3_repo):
    data=p3_repo/'particle_3d/reports/particle3/data'
    audit=json.loads((data/'04_real_contact_velocity_audit.json').read_text())
    assert audit['status']=='PASS' and audit['processed_case_count']==18
    assert audit['no_state_or_velocity_modified'] is True
    assert audit['maximum_residual_to_bound_ratio']<=1
    assert audit['qp_roundoff_factor']==1024 and audit['normal_merge_direction_factor']==256
    assert audit['no_wall_or_velocity_tolerance_change'] is True
    expected=0
    for path in data.glob('*_real_*_states.json'):
        rows=json.loads(path.read_text());expected+=sum(row['contact_state']=='TOUCHING' for row in rows[:-1])
    assert audit['contact_intervals']==expected
