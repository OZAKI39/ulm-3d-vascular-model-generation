"""Public independent WSL finalizer: numerical references, exact geometry, visualization."""
from pathlib import Path
import argparse,sys,json,subprocess,hashlib
ap=argparse.ArgumentParser();ap.add_argument('--root',default=str(Path(__file__).parent));ap.add_argument('--output',default='PASSIVE_MICROBUBBLE_TRANSPORT_V0_VALIDATION.json');ap.add_argument('--resume-verified-local-audits',action='store_true');args=ap.parse_args();root=Path(args.root).resolve();sys.path.insert(0,str(root/'src'))
from finalize_passive_transport_v0 import audit_all,integrity

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def verify_manifest(path):
 for line in Path(path).read_text().splitlines():
  digest,name=line.split('  ',1);assert sha(root/name)==digest,name
assert integrity(root)
if args.resume_verified_local_audits:
 # Postprocessing-only recovery: immutable numerical evidence and already completed
 # independent audits must match the recorded byte identities before any reuse.
 verify_manifest(root/'NUMERICAL_EVIDENCE_SHA256SUMS');verify_manifest(root/'LOCAL_AUDIT_COMPONENT_SHA256SUMS')
 num=json.loads((root/'LOCAL_NUMERICAL_FINALIZER.json').read_text())
else:
 num=audit_all(root);(root/'LOCAL_NUMERICAL_FINALIZER.json').write_text(json.dumps(num,indent=2)+'\n');assert num['status']=='PASS'
 for name in ['local_geometry_audit.py','visualize_passive.py']:
  with (root/(name+'.log')).open('w') as log:p=subprocess.run([sys.executable,'-B',str(root/'scripts'/name),str(root)],stdout=log,stderr=subprocess.STDOUT)
  assert p.returncode==0,name+' failed; inspect log'
geo=json.loads((root/'LOCAL_GEOMETRY_AND_REPLAY_AUDIT.json').read_text());vis=json.loads((root/'VISUALIZATION_PROVENANCE.json').read_text());assert num['status']==geo['status']==vis['status']=='PASS'
for filename,item in vis['artifacts'].items():
 assert sha(root/'visualization'/filename)==item['output_sha256'];assert sha(root/item['plot_script'])==item['plot_script_sha256']
 for source in item['sources']:assert sha(root/source['file'])==source['sha256']
assert integrity(root)
r={'status':'PASS','INDEPENDENT_FINALIZER':'PASS','numerical':num,'exact_geometry_and_every_step_RK2_replay':geo,'VISUALIZATION_GENERATION':vis['status'],'HUMAN_VISUAL_REVIEW':'PENDING','frozen_inputs_unchanged':True,'full_microbubble_transport_model':'PENDING','postprocess_resume_used':args.resume_verified_local_audits,'finalizer_sha256':sha(__file__)};out=Path(args.output);out=out if out.is_absolute() else root/out;out.write_text(json.dumps(r,indent=2)+'\n');print('INDEPENDENT_FINALIZER=PASS',str(out),flush=True)
