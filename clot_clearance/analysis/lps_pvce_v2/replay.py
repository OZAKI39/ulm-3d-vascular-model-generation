"""Read-only compact verification. Never launch any equilibrium or overwrite data."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"verification/lps_pvce_continuum_and_near_onset_audit_v2/20261011T142003Z"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sys.argv[1:]==['--verify-only'],'Only --verify-only is implemented; new solves require a new frozen run.'
m=json.loads((BASE/'FILE_MANIFEST.json').read_text());present=missing=0
for r in m['files']:
 p=ROOT/r['path']
 if not p.exists():
  assert p.name=='body.data','Missing required compact artifact: '+r['path'];missing+=1;continue
 assert sha(p)==r['sha256'],r['path'];present+=1
E=json.loads((BASE/'ENERGY_AND_EQUILIBRIUM_AUDIT.json').read_text())
for r in E['comparison']:
 assert abs((r['U0_J']-r['U1_J'])/3e-8-r['G_two_J_m2'])<1e-19
 assert abs((3*r['U0_J']-4*r['U1_J']+r['U2_J'])/6e-8-r['G_three_J_m2'])<1e-19
s=json.loads((BASE/'A0_FIRST_FAILURE_SEARCH.json').read_text());assert s['critical_opening_bracket_m'] is None and s['first_location']=='UNDETERMINED';assert s['admission']['high_A1_A2']=='NOT_RUN';assert not (BASE/'NEAR_ONSET_PVCE_RESULTS.json').exists()
print(json.dumps(dict(manifest_verified=present,server_only_reconstructible_body_files_missing=missing,energy_differences_recomputed=True,high_PVCE_correctly_absent=True,new_solves=0)))
