import json,subprocess
from particle_3d.audit import sha256
from particle_3d.particle6_checkpoint import DEPENDENCIES

def test_particle5_dependency_locked(repo):
    p=repo/'particle_3d/reports/particle5'
    v=json.loads((p/'PARTICLE5_VALIDATION.json').read_text())
    assert v['git_commit']==DEPENDENCIES[5]
    for name,digest in v['source_sha256'].items():assert sha256(repo/name)==digest,name
    a=json.loads((p/'PARTICLE5_MANUAL_ACCEPTANCE.json').read_text())
    assert a['manual_visual_review']=='PASS' and a['particle6_authorized_to_start']
    for commit in DEPENDENCIES:subprocess.run(['git','merge-base','--is-ancestor',commit,'HEAD'],cwd=repo,check=True)
