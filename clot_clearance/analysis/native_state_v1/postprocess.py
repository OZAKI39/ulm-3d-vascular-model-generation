"""Read-only historical graph reconstruction + five deterministic semantic regressions."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
from pd_clot.lammps_fragmentation.detached_fragment_metrics import graph_components
from pd_clot.lammps_fragmentation.fragment_semantics import component_metrics,classify_detached_v2,LEGACY_FIELDS
from pd_clot.lammps_fragmentation.failure_classifier import classify
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path.cwd();O=R/'verification'/C['task']/C['run']
def tests():
 cases=[('all_attached',[20],[True],0,0),('body_and_singles',[20,1,1,1],[True,False,False,False],3,0),('body_and_34',[50,34],[True,False],34,1),('attached_and_free_tiny',[20,2,1],[True,True,False],1,0),('partial_interface_failure',[20],[True],0,0)]
 reports=[]
 for name,sizes,attached,expected_nodes,ge200 in cases:
  n=sum(sizes);b=[];typ=[];offset=0
  for sz,att in zip(sizes,attached):
   for j in range(offset,offset+sz-1):b.append([j,j+1]);typ.append(1)
   if att:b.append([offset,n]);typ.append(2)
   offset+=sz
  if name=='partial_interface_failure':b.extend([[5,n],[6,n]]);typ.extend([2,2])
  b=np.asarray(b,int);typ=np.asarray(typ);alive=np.ones(len(b),bool)
  if name=='partial_interface_failure':alive[-2:]=False
  lab,count,att=graph_components(n,b,typ,alive);m=component_metrics(count,att,50e-6)
  assert sorted(count.tolist())==sorted(sizes) and count.sum()==n
  assert m['detached_component_count']==sum(not a for a in attached)
  assert abs(m['detached_mass_fraction']-expected_nodes/n)<1e-15
  assert abs(m['detached_mass_kg']-expected_nodes*1.25e-10)<1e-22
  assert m['N_detached_ge200']==ge200 and m['N_detached_ge300']==0
  assert abs(m['attached_mass_kg']+m['detached_mass_kg']-m['total_material_mass_kg'])<1e-22
  if name=='attached_and_free_tiny':assert abs(m['attached_tiny_component_mass_fraction']-2/23)<1e-15
  if name=='partial_interface_failure':assert att.all()
  reports.append({'case':name,'status':'PASS','metrics':m})
 # Explicit 300um threshold and heterogeneous masses: no equal-mass assumption hidden.
 z=component_metrics([113,114],[False,False],50e-6);assert z['N_detached_ge300']==1
 z=component_metrics([2,3],[True,False],50e-6,component_mass=[10.,1.]);assert abs(z['detached_mass_fraction']-1/11)<1e-15
 return reports
p=R/'verification/y_clot_point_cloud_and_early_flow_coupling_v1/cases/C1_GAUSSIAN_A200';geom=R/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/V3_MATERIAL_POINTS_AND_BONDS.npz'
g=np.load(geom);print('geometry',g.files);X=g['X'];b=g['bonds'];typ=g['bond_types'];n=1927
tr=np.load(p/'accepted_trajectory.npz');ev=np.load(p/'broken_bonds.npz');print('trajectory',tr.files,'events',ev.files)
rows=[];alive=np.ones(len(b),bool);cursor=0;events=ev['events'];idx=ev['initial_bond_indices'];labels=[];attachments=[]
oldrows=list(csv.DictReader((p/'fragment_metrics.csv').open()))
for k,step in enumerate(tr['step']):
 while cursor<len(idx) and events[cursor,1]<=step:alive[idx[cursor]]=False;cursor+=1
 lab,count,att=graph_components(n,b,typ,alive);stable=np.full(len(count),n+1);np.minimum.at(stable,lab,np.arange(1,n+1))
 assert np.array_equal(stable[lab],tr['component_id'][k]) and np.array_equal(att[lab],tr['attached'][k])
 legacy=dict(debris_fraction=float(count[count<=2].sum()/n),largest_free_fraction=float(count.max()/n),FI_free=float(1-count.max()/n),resolved_fragment_count=int((count>=3).sum()),median_fragment_nodes=float(np.median(count)),bulk_broken_fraction=float((~alive[typ==1]).mean()),interface_broken_fraction=float((~alive[typ==2]).mean()))
 for key in ['debris_fraction','largest_free_fraction','FI_free','resolved_fragment_count']:
  assert abs(float(oldrows[k][key])-legacy[key])<1e-13
 m=component_metrics(count,att,50e-6);row={'step':int(step),'time_s':float(step*.5e-6),**legacy,**m,'failure_mode_legacy_v1':classify(legacy),'failure_mode_detached_v2':classify_detached_v2({**legacy,**m})}
 row['POSTPROCESS_CLASSIFICATION_SEMANTICS_CHANGED']=row['failure_mode_legacy_v1']!=row['failure_mode_detached_v2'];rows.append(row);labels.append(stable[lab]);attachments.append(att[lab])
with (O/'POSTPROCESS_OLD_NEW_COMPARISON.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
np.savez_compressed(O/'historical_final_graph.npz',position_m=tr['position_m'][-1],X=X,bonds=b,types=typ,alive=alive,material_component=labels[-1],attached=attachments[-1],counts=count,component_attached=att,step=tr['step'][-1])
report={'status':'PASS','tests':tests(),'legacy_fields':LEGACY_FIELDS,'legacy_numeric_values_unchanged_all_frames':True,'new_version':'ATTACHMENT_AWARE_V2','legacy_classifier_retained':'ALL_COMPONENTS_V1','historical_graph_exact_all_frames':True,'frames':len(rows),'new_vs_legacy_label_different_frames':sum(r['POSTPROCESS_CLASSIFICATION_SEMANTICS_CHANGED'] for r in rows),'final':rows[-1],'scope':'New separate classifications are opt-in. Historical CSV/NPZ/AI labels untouched. Generic fragment_analysis ALL metrics intentionally include anchored material, not Y wall anchors. 1–2 node debris unresolved.','source_hashes':{str(q.relative_to(R)):hashlib.sha256(q.read_bytes()).hexdigest() for q in [geom,p/'accepted_trajectory.npz',p/'broken_bonds.npz',p/'fragment_metrics.csv']}}
(O/'FRAGMENT_STATISTICS_SEMANTIC_AUDIT.json').write_text(json.dumps(report,indent=2));print(json.dumps(report['final'],indent=2))
