import json,hashlib,subprocess
def test_locked(repo):
 snapshot=json.loads((repo/'particle_3d/reports/particle6_5/PARTICLE6_5_VALIDATION.json').read_text())
 for name,digest in snapshot['source_sha256'].items():
  assert hashlib.sha256((repo/name).read_bytes()).hexdigest()==digest,name
 subprocess.run(['git','merge-base','--is-ancestor','57dfe7bd15b29aca9a941ab4fd74dae081d9f161','HEAD'],cwd=repo,check=True)
