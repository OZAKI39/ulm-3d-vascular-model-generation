import json,hashlib,subprocess,sys
from pathlib import Path
import numpy as np
from PIL import Image
from particle_3d.particle_shapes import Sphere,EPS
from particle_3d.pair_geometry import pair_gap
from particle_3d.wall_gap import wall_gap
from particle_3d.wall_geometry import WallGeometry
from particle_3d.field import FrozenFEMField


def test_particle3_dependency_and_acceptance_are_locked(p4_repo):
    p3=p4_repo/'particle_3d/reports/particle3'
    v=json.loads((p3/'PARTICLE3_VALIDATION.json').read_text());acceptance=json.loads((p3/'PARTICLE3_MANUAL_ACCEPTANCE.json').read_text())
    assert v['git_commit']=='95e54fe474349c35aaad2e2754aff0bad7108c42'
    assert acceptance['manual_visual_review']=='PASS' and acceptance['particle4_authorized_to_start']
    assert acceptance['real_rbc_passage']=='NOT_ESTABLISHED' and acceptance['physiological_blockage_conclusion']=='NOT_ALLOWED'
    for group in ['source_sha256','data_sha256','figure_sha256']:
        for path,digest in v[group].items():assert hashlib.sha256((p4_repo/path).read_bytes()).hexdigest()==digest,path


def test_scope_has_no_new_physical_contact_parameters(p4_repo):
    scope=json.loads((p4_repo/'particle_3d/reports/particle4/data/00_scope.json').read_text())
    assert scope['collision_induced_deformation'] is False and scope['physical_force_claimed'] is False
    assert scope['particle_particle_lubrication']=='DEFERRED_PARTICLE5'
    assert scope['particle5_started'] is False and scope['no_cfd_executed'] is True
    assert scope['production_particle_timestep_frozen'] is False
    from inspect import signature
    from particle_3d.kinematic_contact import project_contacts
    assert list(signature(project_contacts).parameters)==['shapes','free_velocities','free_omegas','constraints']


def test_real_fem_two_mb_smoke_finite_inside_clear_and_original_sizes(p4_repo):
    data=p4_repo/'particle_3d/reports/particle4/data';rows=json.loads((data/'14_real_two_mb_states.json').read_text())
    summary=json.loads((data/'14_real_two_mb_summary.json').read_text());initial=json.loads((data/'14_initialization.json').read_text())
    sizes=json.loads((data/'14_sonovue_provenance.json').read_text());assert [s['seed'] for s in sizes]==[20260920,20260921]
    assert summary['status']=='PASS' and summary['no_artificial_forces'] is True
    assert all(initial['inside_lumen'].values()) and initial['initial_overlap'] is False
    wall=WallGeometry.from_frozen(p4_repo/'formal_3D_flow_solver/FEM_SimVascular')
    field=FrozenFEMField.from_frozen(p4_repo/'formal_3D_flow_solver/FEM_SimVascular')
    for row in rows:
        for k,p in enumerate(row['particles']):
            assert p['radius_m']==sizes[k]['radius_m']
            assert np.isfinite(p['center_m']+p['velocity_m_s']+p['omega_s_inv']).all()
            assert p['wall_gap_m']>=-p['wall_roundoff_m']
            assert field.sample(p['center_m']).inside_lumen
        assert all(g['gap_m']>=-g['roundoff_budget_m'] for g in row['pair_gaps'])
    for row in [rows[0],rows[len(rows)//2],rows[-1]]:
        shapes=[Sphere(p['center_m'],p['radius_m']) for p in row['particles']]
        g=pair_gap(*shapes);assert abs(g.gap_m-row['pair_gaps'][0]['gap_m'])<=g.roundoff_budget_m
        for shape,p in zip(shapes,row['particles']):
            gap=wall_gap(shape,wall);assert abs(gap.gap_m-p['wall_gap_m'])<=gap.roundoff_m


def test_evidence_intervals_kkt_and_mixed_no_overlap(p4_repo):
    data=p4_repo/'particle_3d/reports/particle4/data'
    for prefix in ['13_mixed_dt0','13_mixed_dt1','13_mixed_dt2','14_real_two_mb']:
        rows=json.loads((data/(prefix+'_states.json')).read_text());ledger=json.loads((data/(prefix+'_intervals.json')).read_text())
        assert abs(sum(r['dt_s'] for r in ledger)-rows[-1]['time_s'])<=16*EPS
        assert all(a['t1_s']==b['t0_s'] for a,b in zip(ledger,ledger[1:]))
        assert all(r['dt_s']>0 and r['accepted'] and r['proof'] for r in ledger)
        for row in rows:
            assert all(g['gap_m']>=-g['roundoff_budget_m'] for g in row['pair_gaps'])
            audit=row['projection']
            if audit:
                assert audit['max_primal_violation_m_s']<=audit['velocity_budget_m_s']
                assert audit['max_dual_violation']<=audit['multiplier_budget']
                assert audit['max_complementarity_residual']<=audit['complementarity_budget']


def test_all_16_figures_and_persisted_source_hashes(p4_repo):
    report=p4_repo/'particle_3d/reports/particle4';figures=sorted((report/'figures').glob('*.png'));assert len(figures)==16
    for i,path in enumerate(figures):
        assert path.name.startswith(f'{i:02d}_')
        with Image.open(path) as im:assert im.width>=1200 and im.height>=650
        source=json.loads((report/f'data/{i:02d}_figure_sources.json').read_text())
        for name,digest in source['source_sha256'].items():assert hashlib.sha256((report/'data'/name).read_bytes()).hexdigest()==digest
        notes=(report/f'{i:02d}_step_notes.md').read_text()
        for label in ['应该看什么','实际看到什么','有没有异常']:assert notes.count(label)==1


def test_plot_reproduction(p4_repo,tmp_path):
    script=p4_repo/'particle_3d/scripts/generate_particle4_report.py'
    run=subprocess.run([sys.executable,'-B',str(script),'--only','2','--output',str(tmp_path)],capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    generated=next(tmp_path.glob('*.png'));original=p4_repo/'particle_3d/reports/particle4/figures'/generated.name
    assert hashlib.sha256(generated.read_bytes()).digest()==hashlib.sha256(original.read_bytes()).digest()


def test_review_retains_manual_pending_and_particle3_limitation(p4_repo):
    report=p4_repo/'particle_3d/reports/particle4';v=json.loads((report/'PARTICLE4_VALIDATION.json').read_text());text=(report/'PARTICLE4_REVIEW.md').read_text()
    assert v['manual_visual_review']=='PENDING_USER_REVIEW' and v['particle5_started'] is False
    assert v['particle3_real_rbc_passage']=='NOT_ESTABLISHED' and v['physical_force_claimed'] is False
    assert sum(line.startswith('## ') for line in text.splitlines())==11
    assert all('figures/'+p.name in text for p in (report/'figures').glob('*.png'))
