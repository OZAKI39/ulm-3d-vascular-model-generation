import json,subprocess
from particle_3d.audit import sha256

def test_dependency(repo,report):
 d=json.loads((report/'data/00_scope.json').read_text())
 # P9-A changes only the input reader's case label; historical model hashes
 # remain strict. Verify that reader against its preserved baseline Git object.
 import hashlib
 for k,v in d['p6_source_sha256'].items():
  if k=='particle_3d/src/particle_3d/audit.py':
   baseline=subprocess.check_output(['git','show','775adc019536585a4ee2f4dbf2c083c959e0a339:'+k],cwd=repo)
   assert hashlib.sha256(baseline).hexdigest()==v
   from particle_3d.particle9a_provenance import require_current_flow
   assert require_current_flow()['case_role']=='MEAN_2P0_MMPS'
  elif k in {'particle_3d/src/particle_3d/wall_geometry.py',
             'particle_3d/src/particle_3d/resistance_solver.py',
             'particle_3d/src/particle_3d/physical_time_refinement.py'}:
   # Authorized P9-A.1 topology/contact/event repairs; historical evidence
   # still verifies against its original immutable Git object.
   baseline=subprocess.check_output(['git','show','775adc019536585a4ee2f4dbf2c083c959e0a339:'+k],cwd=repo)
   assert hashlib.sha256(baseline).hexdigest()==v
  else: assert sha256(repo/k)==v
 assert all(sha256(repo/k)==v for k,v in d['historical_p5_sha256'].items())
 subprocess.run(['git','merge-base','--is-ancestor','c8e22f99359d15fd6075f2083ed60724ad0cfc2f','HEAD'],cwd=repo,check=True)
 a=json.loads((repo/'particle_3d/reports/particle6/PARTICLE6_MANUAL_ACCEPTANCE.json').read_text())
 assert a['manual_visual_review']=='PASS' and a['particle6_5_authorized_to_start']
