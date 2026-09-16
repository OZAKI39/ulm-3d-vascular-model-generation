from pathlib import Path
import json,hashlib,csv
import numpy as np,h5py
from scipy.spatial import cKDTree
import vtk
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();S=json.loads((R/'validation/INDEPENDENT_FINALIZER.json').read_text());V=json.loads((R/'VISUALIZATION_PROVENANCE.json').read_text())
for section in ['data_sha256','PNG_sha256','VTP_sha256']:
 for name,h in V[section].items():assert sha(R/name)==h,(section,name)
assert sha(R/'scripts/visualize_fixed_multiblob_wall.py')==V['source_script_sha256'];assert len(V['PNG_sha256'])==12 and len(V['VTP_sha256'])==4
assert V['HUMAN_VISUAL_REVIEW']=='PENDING'
for n in ['PLANAR_CONVERGENCE','LATERAL_PHASE_SENSITIVITY','ANISOTROPY','LEAKAGE_SCREENING','MATRIX_PROPERTIES','CYLINDER_CURVATURE','REAL_STL_PATCH_AUDIT','PERFORMANCE_SCALING']:
 with (R/(n+'.csv')).open() as f:assert len(list(csv.DictReader(f)))>0,n
geometries=0;queries=0;failed=0;minimum_gap=float('inf')
with h5py.File(R/'FIXED_MULTIBLOB_WALL_AUDIT.h5') as f:
 for name,root in f.items():
  if not name.startswith('RAW_'):continue
  for gn,g in root.items():
   x=g['coords'][:];a=g['radii'][:];assert len(x)==g.attrs['N_wall'] and np.all(a==g.attrs['beta'])
   if len(x)>1:
    minimum=float((cKDTree(x).query(x,k=2)[0][:,1]-2*a).min());assert minimum>=-1e-12,(name,gn,minimum);minimum_gap=min(minimum_gap,minimum)
   geometries+=1
   for qn,q in g.items():
    if not isinstance(q,h5py.Group):continue
    t=q['targets'][:];ar=q['target_radii'][:];assert np.min(np.linalg.norm(t[:,None]-x[None],axis=2)-ar[:,None]-a[None])>=-1e-12
    if q.attrs['status']!='COMPLETED':failed+=1;continue
    assert np.isfinite(q['R_scaled'][:]).all()
    if q['R_scaled'].shape==(6,6):
     queries+=1;assert q['force_responses_256'].shape==(256,6) and q['power_10000'].shape==(10000,)
assert queries==S['completed_response_queries'] and failed==S['failed_queries']
vtpcounts={}
for p in R.glob('*.vtp'):
 rd=vtk.vtkXMLPolyDataReader();rd.SetFileName(str(p));rd.Update();d=rd.GetOutput();vtpcounts[p.name]=d.GetNumberOfPoints();assert d.GetFieldData().GetAbstractArray('STATUS') is not None
 assert all(d.GetPointData().GetArray(k) is not None for k in ['object_type','blob_radius','layer','patch_id','bubble_id','gap','matrix_error','validity_flags'])
 if p.name!='PLANAR_WALL_BLOBS.vtp':assert d.GetNumberOfPoints()==0
assert vtpcounts['PLANAR_WALL_BLOBS.vtp']>0
res=dict(status='PASS',raw_geometry_groups=geometries,completed_6x6_queries=queries,failed_queries_preserved=failed,minimum_wall_wall_gap_over_a=minimum_gap,PNG_count=12,VTP_counts=vtpcounts,visualization_hashes_match=True,human_visual_review='PENDING',source_and_reference_identity=json.loads((R/'validation/IMMUTABILITY_AFTER.json').read_text())['status'],download_identity=json.loads((R/'validation/REMOTE_TO_WSL_INTEGRITY.json').read_text())['status']);(R/'validation/FINAL_ARTIFACT_CHECK.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
