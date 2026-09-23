"""Post-run checks and structured extrema; executed on the same audit server."""
from pathlib import Path
import sys,json,hashlib,xml.etree.ElementTree as ET
import numpy as np
from run_audit import dump,sha

def main(config):
 cfg=json.loads(Path(config).read_text());out=Path(cfg['output']);base=Path(cfg['input_base'])
 inv=json.loads(Path(cfg['inventory']).read_text());receipts=json.loads((out/'trajectory_audit_receipts.json').read_text())
 extrema={};dataset_counts={};state_count=0;initial=0;maxgap=0.;outputs={}
 for r in receipts:
  state_count+=r['count']
  if not r['count']:continue
  maxgap=max(maxgap,r['max_gap_error_m']);f=Path(r['file'])
  with np.load(f) as z:
   initial+=int(np.sum(z['initial_diagnostic']));ds=r['dataset'];counts=dataset_counts.setdefault(ds,{'total':0,'far':0,'near':0,'handoff':0})
   counts['total']+=r['count']
   for name,key in [('far','region_far_field'),('near','region_near_wall'),('handoff','region_handoff')]:counts[name]+=int(np.sum(z[key]))
   for metric in ['CANDIDATE_SAFFMAN_MAGNITUDE','slip_speed_m_s','shear_s_inv','lift_drag_ratio','lift_lubrication_ratio']:
    a=z[metric];ok=np.flatnonzero(np.isfinite(a))
    if not len(ok):continue
    k=int(ok[np.argmax(a[ok])]);v=float(a[k])
    if metric not in extrema or extrema[metric]['value']<v:
     extrema[metric]=dict(dataset=ds,id=r['id'],index=k,value=v,x_m=z['evaluation_x_m'][k].tolist(),time_s=float(z['evaluation_time_s'][k]),
       radius_m=float(z['radius_m'][k]),gap_m=float(z['h_m'][k]),h_over_a=float(z['h_over_a'][k]),state=int(z['state'][k]),
       shear_s_inv=float(z['shear_s_inv'][k]),slip_m_s=float(z['slip_speed_m_s'][k]),wall_margin=float(z['wall_margin'][k]),
       candidate_N=float(z['CANDIDATE_SAFFMAN_MAGNITUDE'][k]),drag_N=float(z['drag_N'][k]),lubrication_N=float(z['lubrication_N'][k]),
       regions=[key[7:] for key in z.files if key.startswith('region_') and z[key][k]])
 assert state_count==sum(r['sample_count'] for r in inv['records'] if not r['duplicate_of'])
 changed=[]
 for r in inv['records']:
  if r['duplicate_of']:continue
  p=base/r['remote_path']
  if sha(p)!=r['sha256'] or sha(p.with_suffix('.json'))!=r['metadata_sha256']:changed.append(str(p))
 assert not changed
 xml=ET.parse(base/'new_solver.xml');mu=float(xml.findtext('.//Viscosity/Value'));rho=float(xml.findtext('.//Add_equation/Density'))
 assert mu==.00345312 and rho==1056
 dump(out/'post_verification.json',dict(status='PASS',compute_host=__import__('socket').gethostname(),state_count=state_count,
   initial_diagnostic_rows=initial,accepted_interval_rows=state_count-initial,maximum_gap_reconstruction_error_m=maxgap,
   dataset_state_counts=dataset_counts,extrema=extrema,changed_inputs=changed,new_flow_density=rho,new_flow_viscosity=mu))
 for p in out.rglob('*'):
  if p.is_file() and p.name!='OUTPUT_SHA256.json':outputs[str(p.relative_to(out))]=sha(p)
 dump(out/'OUTPUT_SHA256.json',outputs)
 print('POST_VERIFICATION_PASS',state_count,len(outputs))
if __name__=='__main__':main(sys.argv[1])
