from pathlib import Path
import importlib.util,json,subprocess,sys
import numpy as np
import pytest
from PIL import Image
from particle_3d.audit import sha256
from particle_3d.particle_shapes import Sphere,EPS
from particle_3d.particle5_cases import P4_COMMIT


def load(repo,name):return json.loads((repo/'particle_3d/reports/particle5/data'/f'{name}.json').read_text())


def test_particle4_dependency_locked(repo):
    validation=json.loads((repo/'particle_3d/reports/particle4/PARTICLE4_VALIDATION.json').read_text())
    assert validation['git_commit']==P4_COMMIT
    for name,digest in validation['source_sha256'].items():assert sha256(repo/name)==digest,name
    acceptance=json.loads((repo/'particle_3d/reports/particle4/PARTICLE4_MANUAL_ACCEPTANCE.json').read_text())
    assert acceptance['manual_visual_review']=='PASS' and acceptance['particle5_authorized_to_start']
    subprocess.run(['git','merge-base','--is-ancestor',P4_COMMIT,'HEAD'],cwd=repo,check=True)


def test_particle5_scope(repo):
    scope=load(repo,'00_scope')
    assert scope['resistance_formulation']=='AFFINE_BLOCK_RESISTANCE_BALANCE_V0'
    assert scope['matrix_dimension']=='6N' and scope['nonspherical_lubrication_model']=='NOT_FROZEN'
    assert not scope['production_lubrication_cutoff_frozen'] and not scope['production_particle_timestep_frozen']
    assert scope['production_wall_lubrication_cutoff'] is scope['production_pair_lubrication_cutoff'] is None
    assert scope['validation_nearfield_ratio_max']==.01
    assert scope['real_rbc_hydrodynamics']=='DEFERRED' and scope['particle3_real_rbc_passage']=='NOT_ESTABLISHED'
    assert scope['no_cfd_executed'] and not scope['particle6_started']


def test_real_near_wall_mb_lubrication(repo):
    from particle_3d.field import FrozenFEMField
    from particle_3d.wall_geometry import WallGeometry
    from particle_3d.wall_gap import wall_gap
    fem=repo/'formal_3D_flow_solver/FEM_SimVascular';field=FrozenFEMField.from_frozen(fem);wall=WallGeometry.from_frozen(fem)
    replays=load(repo,'09_real_near_wall')
    assert replays[0]['free_normal_velocity_m_s']>0 and replays[1]['free_normal_velocity_m_s']<0
    for row in replays:
        source=repo/row['source_path'];assert sha256(source)==row['source_sha256']
        state=json.loads(source.read_text())[row['source_row_index']]
        p=next(p for p in state['particles'] if p['particle_id']==row['original_particle']['particle_id'])
        assert row['center_m']==p['center_m'] and row['radius_m']==p['radius_m']
        shape=Sphere(p['center_m'],p['radius_m']);g=wall_gap(shape,wall);sample=field.sample(shape.center_m)
        assert g.gap_m==row['gap_m'] and g.gap_m>g.roundoff_m and row['gap_ratio']<=.01
        np.testing.assert_array_equal(row['free'][:3],sample.velocity_m_s)
        assert sample.inside_lumen
        vn=np.asarray(sample.velocity_m_s)@g.normal_inward;expected=vn*g.gap_m/(g.gap_m+shape.radius_m)
        budget=64*EPS*row['solver']['condition_estimate']*np.linalg.norm(sample.velocity_m_s)
        assert abs(row['lubricated_normal_velocity_m_s']-expected)<budget
        assert 0<row['attenuation_ratio']<.01 and row['tangential_difference_m_s']<budget


def test_real_two_mb_resistance_smoke(repo):
    states=load(repo,'10_real_two_mb_states');summary=load(repo,'10_real_summary');initial=load(repo,'10_initialization')
    p4=json.loads((repo/initial['source']).read_text());assert initial['source_sha256']==sha256(repo/initial['source'])
    assert states[0]==p4['initial_state']
    assert summary['requested_dt_s']==p4['validation_dt_s'] and summary['final_time_s']==p4['horizon_s']
    assert summary['pair_lubrication_status']=='PAIR_LUBRICATION_NOT_ACTIVE_IN_THIS_REAL_SMOKE'
    assert summary['minimum_pair_gap_ratio']>.01 and summary['pair_active_block_count']==0
    ledger=load(repo,'10_real_two_mb_intervals');assert len(ledger)==len(states)-1
    for index,state in enumerate(states):
        for p,original in zip(state['particles'],states[0]['particles']):
            assert p['radius_m']==original['radius_m']
            assert p['wall_gap_m']>=-p['wall_roundoff_m']
            assert np.isfinite(p['center_m']+p['velocity_m_s']+p['omega_s_inv']).all()
        assert all(g['gap_m']>=-g['roundoff_budget_m'] for g in state['pair_gaps'])
        if index:
            previous=states[index-1];dt=state['time_s']-previous['time_s']
            assert dt>0 and abs(dt-ledger[index-1]['dt_s'])<1e-17
            for old,new in zip(previous['particles'],state['particles']):
                expected=np.array(old['center_m'])+dt*np.array(new['velocity_m_s'])
                np.testing.assert_allclose(new['center_m'],expected,rtol=512*EPS,atol=0)
            assert state['projection']['dissipation']['total']>=0
            for b in state['projection']['pair_eligibility']:
                assert not b['active'] and b['reason']=='OUTSIDE_VALIDATION_ASYMPTOTIC_RANGE'
    assert abs(sum(l['dt_s'] for l in ledger)-summary['final_time_s'])<1e-15


def test_timestep_independent_implicit_ode_convergence(repo):
    rows=load(repo,'11_timestep')
    for kind in ['wall','pair']:
        errors=[abs(r['implicit_ode_defect_m']) for r in rows if r['kind']==kind]
        assert errors[2]<errors[1]<errors[0]
        assert 1.8<errors[0]/errors[1]<2.2 and 1.8<errors[1]/errors[2]<2.2


def test_particle5_plot_generation(repo,tmp_path):
    script=repo/'particle_3d/scripts/generate_particle5_report.py'
    subprocess.run([sys.executable,'-B',str(script),'--output',str(tmp_path),'--no-notes'],cwd=repo,check=True)
    manifest=json.loads((repo/'particle_3d/reports/particle5/FIGURE_MANIFEST.json').read_text())
    assert len(manifest)==12
    for row in manifest:
        assert sha256(tmp_path/row['figure'])==row['sha256']
        for source,digest in row['sources'].items():assert sha256(repo/source)==digest
        with Image.open(tmp_path/row['figure']) as image:assert image.width>=1500 and image.height>=700


def test_particle5_review_artifacts(repo):
    report=repo/'particle_3d/reports/particle5'
    review=(report/'PARTICLE5_REVIEW.md').read_text()
    assert sum(line.startswith('## ') for line in review.splitlines())==12
    for index in range(12):
        note=(report/f'{index:02d}_step_notes.md').read_text()
        for text in ['应该看什么','实际看到什么','有没有异常']:assert text in note and text in review
    v=json.loads((report/'PARTICLE5_VALIDATION.json').read_text())
    for key in ['max_stokes_resistance_error','max_wall_lubrication_scaling_error','max_pair_lubrication_scaling_error',
                'max_analytic_velocity_error','matrix_symmetry_error','minimum_dissipation','resistance_solver_residual','max_condition_number']:
        assert np.isfinite(v[key])
    assert v['manual_visual_review']=='PENDING_USER_REVIEW' and not v['particle6_started']
    assert v['sparse_all_pairs_equivalence'] and v['real_rbc_hydrodynamics']=='DEFERRED'
    for key in ['source_sha256','data_sha256','figure_sha256']:
        for name,digest in v[key].items():assert sha256(repo/name)==digest
