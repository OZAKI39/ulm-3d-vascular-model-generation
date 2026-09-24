import json,subprocess
from particle_3d.audit import sha256

def test_dependency(repo,report):
 d=json.loads((report/'data/00_scope.json').read_text())
 assert all(sha256(repo/k)==v for k,v in d['p6_source_sha256'].items())
 assert all(sha256(repo/k)==v for k,v in d['historical_p5_sha256'].items())
 subprocess.run(['git','merge-base','--is-ancestor','c8e22f99359d15fd6075f2083ed60724ad0cfc2f','HEAD'],cwd=repo,check=True)
 a=json.loads((repo/'particle_3d/reports/particle6/PARTICLE6_MANUAL_ACCEPTANCE.json').read_text())
 assert a['manual_visual_review']=='PASS' and a['particle6_5_authorized_to_start']
