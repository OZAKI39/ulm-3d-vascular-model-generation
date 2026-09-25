"""Reevaluate numerical regression evidence; preserve required upstream-source gate failure."""
from pathlib import Path
import argparse,json,subprocess,sys,hashlib,os

def main(root,out):
 root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 results={}
 for stage,script in [('engine','finalize_lammps_particle_engine.py'),('coupling','finalize_palabos_lammps_coupling.py'),('passive','finalize_passive_transport_v0.py')]:
  r=root/'migration'/stage;target=out/(stage+'_INDEPENDENT.json')
  cmd=[sys.executable,'-B',str(r/'src'/script),'--root',str(r),'--output',str(target)]
  p=subprocess.run(cmd,env=env,capture_output=True,text=True)
  if not target.exists():raise RuntimeError(p.stdout+p.stderr)
  data=json.loads(target.read_text());results[stage]={'status':data['status'],'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr,'output_sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
 from finalize_source_correction_audit import audit
 correction=audit(root);(out/'SOURCE_CORRECTION_RECOMPUTED.json').write_text(json.dumps(correction,indent=2)+'\n')
 required=json.loads((root/'contracts/LAMMPS_2026_MIGRATION_CONTRACT.json').read_text())
 source_ok=correction['status']=='PASS'
 result={'status':'PASS' if source_ok and all(v['status']=='PASS' for v in results.values()) else 'FAIL','regressions':results,'source_correction_gate':correction['status'],'target_commit':required['target']['target_commit'],'promotion_allowed':source_ok and all(v['status']=='PASS' for v in results.values()),'phase_C_allowed':False,'scope':'Independent numeric results + upstream correction gate. Phase B integrity requirements additionally apply before any Phase C.'}
 (out/'INDEPENDENT_MIGRATION_EVALUATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
 # A faithful diagnosis of a failed source gate is a successfully executed audit, not scientific PASS.
 return 0
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--output-dir',required=True);a=p.parse_args();raise SystemExit(main(a.root,a.output_dir))
