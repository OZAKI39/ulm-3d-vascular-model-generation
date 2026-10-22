import importlib.util,json,subprocess,sys
from pathlib import Path
from PIL import Image
from particle_3d.audit import sha256


def test_particle3_all_13_figures_and_source_records(p3_repo):
    report=p3_repo/'particle_3d/reports/particle3';figs=sorted((report/'figures').glob('*.png'))
    assert len(figs)==13
    for i,path in enumerate(figs):
        assert path.name.startswith(f'{i:02d}_')
        with Image.open(path) as img:assert img.width>=1200 and img.height>=650
        source=json.loads((report/f'data/{i:02d}_figure_sources.json').read_text())
        assert source['figure']==path.name and source['manual_visual_review']=='PENDING_USER_REVIEW'
        for relative,digest in source['source_sha256'].items():assert sha256(report/'data'/relative)==digest
        notes=(report/f'{i:02d}_step_notes.md').read_text()
        for term in ['应该看什么','实际看到什么','有没有异常']:assert notes.count(term)==1


def test_particle3_plot_generation_is_reproducible(p3_repo,tmp_path):
    command=[sys.executable,'-B',str(p3_repo/'particle_3d/scripts/generate_particle3_report.py'),'--only','2','--output',str(tmp_path)]
    run=subprocess.run(command,cwd=p3_repo,capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    generated=tmp_path/'02_sphere_wall_gap_validation.png';saved=p3_repo/'particle_3d/reports/particle3/figures'/generated.name
    assert sha256(generated)==sha256(saved)


def test_particle3_review_preserves_user_review_and_scope(p3_repo):
    report=p3_repo/'particle_3d/reports/particle3';v=json.loads((report/'PARTICLE3_VALIDATION.json').read_text())
    assert v['stage']=='Particle-3'
    assert v['manual_visual_review']=='PENDING_USER_REVIEW'
    assert v['production_wall_model_status']=='V0_VALIDATION_ONLY_PENDING_USER_REVIEW'
    assert v['distribution_modified'] is False and v['no_cfd_executed'] is True
    assert v['particle4_started'] is False and v['production_particle_timestep_frozen'] is False
    assert v['wall_lubrication']=='DEFERRED_PARTICLE5' and v['transit_time_penalty']=='DEFERRED_PARTICLE5'
    text=(report/'PARTICLE3_REVIEW.md').read_text()
    assert sum(line.startswith('## ') for line in text.splitlines())==12
    assert all(f'figures/{p.name}' in text for p in sorted((report/'figures').glob('*.png')))
