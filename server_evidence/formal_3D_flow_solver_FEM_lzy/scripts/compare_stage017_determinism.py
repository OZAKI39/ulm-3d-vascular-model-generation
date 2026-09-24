#!/usr/bin/env python3
"""Compare independent executions, ignoring IDs but never geometric differences."""
import json,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from fem3d.audit import write_json,timestamp
out=ROOT/'outputs/stage01_7';read=lambda p:json.loads(p.read_text())
a=read(out/'optimizer_result.json');b=read(out/'determinism/optimizer_result.json');checks={}
def canonical(d,key):
 p=d['points_m'];cells=d[key]
 # Only used vertices matter for individual port trial artifacts.
 ids=np.unique(cells);points=p[ids];order=np.lexsort(points[:,::-1].T);sorted_points=points[order]
 remap=np.full(len(p),-1,dtype=np.int64);remap[ids[order]]=np.arange(len(ids))
 conn=np.sort(remap[cells],axis=1);conn=conn[np.lexsort(conn[:,::-1].T)]
 h=hashlib.sha256();h.update(sorted_points.astype('<f8').tobytes());h.update(conn.astype('<i8').tobytes());return h.hexdigest()
for name,search in a['initial_surface_searches'].items():
 other=b['initial_surface_searches'][name]
 checks[name+'_selected_trial']=search['selected_trial']==other['selected_trial']
 checks[name+'_H_values']=[t['H'] for t in search['trials']]==[t['H'] for t in other['trials']]
 for t,u in zip(search['trials'],other['trials']):
  x=np.load(out/t['path']/'cap_mesh.npz');y=np.load(out/'determinism'/u['path']/'cap_mesh.npz')
  checks[name+'_'+t['trial']+'_geometry_connectivity']=canonical(x,'cap_triangles')==canonical(y,'cap_triangles')
checks['iteration_count']=len(a['iterations'])==len(b['iterations'])
for x,y in zip(a['iterations'],b['iterations']):
 name=x['iteration'];r=np.load(out/name/'mesh/volume_mesh.npz');s=np.load(out/'determinism'/name/'mesh/volume_mesh.npz')
 checks[name+'_all_tetra_canonical']=canonical(r,'tetra')==canonical(s,'tetra')
 checks[name+'_volume_QC']=x['volume']['quality']==y['volume']['quality']
 checks[name+'_proxy']=x['volume']['proxy']==y['volume']['proxy']
 checks[name+'_decision']=x['decision']==y['decision']
checks['selected_iteration']=a['selected_iteration']==b['selected_iteration'];checks['termination_reason']=a['termination_reason']==b['termination_reason']
result={'status':'PASS' if all(checks.values()) else 'FAIL','timestamp':timestamp(),'checks':checks,'method':'Two independent Gmsh optimizer executions. Lexicographically sort used physical vertex coordinates, renumber and sort each element connectivity, SHA256 canonical coordinate/connectivity bytes; H values and all volume quality diagnostics compared exactly.','new_volume_meshes_total':len(read(out/'volume_mesh_ledger.json')['entries'])}
write_json(out/'determinism_comparison.json',result);print(json.dumps(result,indent=2));assert result['status']=='PASS'
