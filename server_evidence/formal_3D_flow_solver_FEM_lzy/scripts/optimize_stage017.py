#!/usr/bin/env python3
"""Bounded adaptive port optimizer; all source/config owned by WSL."""
import argparse,json,sys,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from fem3d.audit import sha256,write_json,timestamp
from fem3d.adaptive_port import surface_search,volume_acceptance,feedback,select_iteration
from fem3d.adaptive_surface import vascular_trial,combine_ports
from fem3d.adaptive_volume import generate_volume,dolfinx_loadability
from fem3d.adaptive_qc import measure_volume

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run',choices=['production','determinism'],default='production');a=ap.parse_args()
 root=ROOT/'outputs/stage01_7';out=root if a.run=='production' else root/'determinism';out.mkdir(parents=True,exist_ok=True)
 if (out/'optimizer_result.json').exists():raise RuntimeError('Completed optimizer run already exists')
 inp=ROOT/'inputs/stage01_7';read=lambda p:json.loads(p.read_text())
 policy=read(inp/'acceptance_policy.json');lock=read(inp/'freeze_lock.json');baseline=read(inp/'baseline_recomputed.json');contract=read(inp/'planar_port_contract_v2.json')
 for filename,key in [('acceptance_policy.json','policy_sha256'),('planar_port_contract_v2.json','contract_sha256'),('baseline_recomputed.json','baseline_recomputed_sha256')]:assert sha256(inp/filename)==lock[key]
 assert sha256(ROOT/'inputs/stage01/tagged_surface_si.npz')==lock['source_surface_sha256']
 assert sha256(ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz')==lock['baseline_mesh_sha256']
 assert read(root/'synthetic/results.json')['status']=='PASS','Synthetic tests must precede vascular optimization'
 source=dict(np.load(ROOT/'inputs/stage01/tagged_surface_si.npz'));searches={};selected={};histories={};records=[];termination=None
 def save_trial(name,H,index,reason):
  path=out/'surface_trials'/name/f'trial_{index:02d}'
  if path.exists():raise RuntimeError('Trial exists; no overwrite or hidden retry')
  path.mkdir(parents=True);data,r=vascular_trial(source,contract['ports'][name],H,policy)
  np.savez_compressed(path/'cap_mesh.npz',**data);r.update(trial=f'trial_{index:02d}',path=str(path.relative_to(out)),timestamp=timestamp(),reason=reason,policy_sha256=lock['policy_sha256'])
  r['artifact_sha256']=sha256(path/'cap_mesh.npz');write_json(path/'trial.json',r)
  print(name,r['trial'],'H',H,'triangles',r['triangle_count'],'geometry',r['geometry_status'],'quality',r['quality']['status'] if r['quality'] else 'NOT EVALUATED',flush=True)
  return r
 for name,port in contract['ports'].items():
  edges=source['points_m'][np.asarray(port['rim_edge_ids'])];hrim=float(np.median(np.linalg.norm(edges[:,1]-edges[:,0],axis=1)));req=float(np.sqrt(port['formal_projected_area_m2']/np.pi))
  result=surface_search(lambda H,index,reason:save_trial(name,H,index,reason),baseline['h_volume_m'],hrim,req,policy)
  searches[name]=result;histories[name]=list(result['trials']);write_json(out/'surface_search.json',{'timestamp':timestamp(),'ports':searches})
  if result['status']!='PASS':termination='SURFACE_OPTIMIZER_FAILED';break
  selected[name]=next(t for t in result['trials'] if t['trial']==result['selected_trial'])
 changed=None;reason='Initial combination of independent coarsest passing port surfaces';before={}
 limit=policy['limits']['maximum_volume_iterations']
 for number in range(limit):
  if termination:break
  name=f'iteration_{number:02d}';base=out/name
  # Shared ledger counts attempts, even if a mesher subsequently fails.
  ledger=root/'volume_mesh_ledger.json';entries=read(ledger)['entries'] if ledger.exists() else []
  if len(entries)>=policy['limits']['maximum_total_volume_meshes']:
   termination='MAX_TOTAL_VOLUME_MESHES_REACHED';break
  data={p:dict(np.load(out/t['path']/'cap_mesh.npz')) for p,t in selected.items()}
  surface,sq=combine_ports(source,contract,data,policy)
  (base/'surface').mkdir(parents=True,exist_ok=True);np.savez_compressed(base/'surface/tagged_surface_si.npz',**surface);write_json(base/'qc/surface_invariants.json',sq)
  if sq['status']!='PASS':termination='COMBINED_SURFACE_FAILED';break
  current={p:t['H'] for p,t in selected.items()}
  row={'iteration':name,'changed_port':changed,'reason':reason,'parameters_before':before,'parameters_after':current,'surface_status':sq['status'],'surface':sq,
       'per_port':{p:{'H':t['H'],'h_rim':t['h_rim'],'R_eq':t['R_eq'],'trial':t['trial'],'triangle_count':t['triangle_count'],'q_min':t['quality']['q_tri']['minimum'],'P5':t['quality']['q_tri']['P5'],'median':t['quality']['q_tri']['median'],'trial_artifact_sha256':t['artifact_sha256']} for p,t in selected.items()}}
  write_json(base/'metadata/iteration_inputs.json',row)
  entries.append({'run':a.run,'iteration':name,'timestamp':timestamp(),'path':str(base.relative_to(root))});write_json(ledger,{'maximum':policy['limits']['maximum_total_volume_meshes'],'entries':entries})
  try:
   volume=generate_volume(surface,base,policy);qc=measure_volume(volume,surface,source,contract,policy)
   load=dolfinx_loadability(base/'mesh/fluid.msh');write_json(base/'qc/dolfinx_loadable.json',load)
   acc=volume_acceptance(qc,baseline,policy,load['status']=='PASS');write_json(base/'qc/volume_quality.json',qc);write_json(base/'qc/acceptance.json',acc)
  except Exception as exc:
   row.update(error=repr(exc),decision={'decision':'STOP_INVALID_VOLUME','termination_reason':'VOLUME_VALIDITY_FAILED','stop':True},acceptance={'status':'FAIL'});records.append(row);termination='VOLUME_VALIDITY_FAILED';write_json(base/'qc/failed_iteration.json',row);break
  decision=feedback(number,qc,acc,policy);row.update(volume=qc,acceptance=acc,decision=decision);records.append(row)
  write_json(base/'qc/iteration_result.json',row);write_json(out/'adaptive_decision_log.json',{'policy_sha256':lock['policy_sha256'],'iterations':records})
  print(name,'tetra',qc['proxy']['N_tetra'],'cost',acc['C_P2'],acc['C_tetra'],'low',qc['quality']['total_below_0_1'],qc['quality']['cap_adjacent_below_0_1'],decision,flush=True)
  if decision['stop']:termination=decision['termination_reason'];break
  changed=decision['changed_port'];before=current;reason=decision['reason']
  if len(histories[changed])>=policy['limits']['maximum_surface_trials_per_port']:
   termination='MAX_SURFACE_TRIALS_REACHED';row['decision']={'decision':'STOP_MAX_SURFACE_TRIALS','termination_reason':termination,'stop':True,'requested_port':changed};write_json(base/'qc/iteration_result.json',row);break
  trial=save_trial(changed,current[changed]/decision['H_divisor'],len(histories[changed]),'VOLUME_FEEDBACK_REFINE_ONLY_'+changed.upper());histories[changed].append(trial)
  if trial['geometry_status']!='PASS' or trial['quality']['status']!='PASS':termination='FEEDBACK_SURFACE_FAILED';break
  selected[changed]=trial
 if termination is None:termination='MAX_ADAPTIVE_ITERATIONS_REACHED'
 chosen=select_iteration(records) if all('volume' in r for r in records) else None
 result={'timestamp':timestamp(),'status':'FEASIBLE_PENDING_ROUNDTRIP' if chosen else 'FAIL','run':a.run,'selected_iteration':chosen,'termination_reason':termination,
  'initial_surface_searches':searches,'port_trial_counts':{p:len(v) for p,v in histories.items()},'iterations':records,'policy_sha256':lock['policy_sha256'],'contract_sha256':lock['contract_sha256'],'gpu_used':False,'fem_solved':False,'stage3_started':False}
 write_json(out/'adaptive_decision_log.json',result);write_json(out/'optimizer_result.json',result)
 if chosen and a.run=='production':
  dest=out/'selected';shutil.copytree(out/chosen,dest);shutil.copyfile(inp/'planar_port_contract_v2.json',dest/'planar_port_contract_v2.json')
  write_json(dest/'metadata/selection.json',{'selected_iteration':chosen,'reason':termination,'policy_sha256':lock['policy_sha256'],'contract_sha256':lock['contract_sha256']})
 print('Optimizer stopped:',termination,'selected:',chosen,flush=True)
if __name__=='__main__':main()
