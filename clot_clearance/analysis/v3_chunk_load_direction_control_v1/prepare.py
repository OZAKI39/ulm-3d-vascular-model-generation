from common import *
import shutil,time,pyvista as pv,vtk
from scipy.spatial import cKDTree
from types import SimpleNamespace
def main():
 start=time.perf_counter();assert not (O/'PREPARED.json').exists()
 for p in [O,F,O/'provenance',O/'full',O/'preflight',O/'source_identity']:p.mkdir(parents=True,exist_ok=True)
 manifest=json.loads((B/'PAYLOAD_MANIFEST.json').read_text())
 for name,item in manifest['files'].items():assert sha(R/name)==item['sha256'],name
 prev=json.loads((PREV/'PAYLOAD_MANIFEST.json').read_text())
 for name,item in prev['files'].items():assert sha(R/name)==item['sha256'],name
 freeze=json.loads((B/'V3_PILOT_FREEZE.json').read_text())
 for name,h in freeze['driver_sha256'].items():assert sha(HIST_A/name)==h
 for name,h in freeze['input_sha256'].items():
  assert sha(B/'provenance'/name)==h;shutil.copy2(B/'provenance'/name,O/'provenance'/name)
 assert sha(LIB)==freeze['library_sha256'];ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);X=ref['X'];ids=np.flatnonzero(ref['areas']>0);assert n==1927 and len(X)==3488 and len(ids)==453 and np.all(ids<n)
 geom=load_npz(R/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/V3_MATERIAL_POINTS_AND_BONDS.npz')
 for k in ['X','bonds','bond_types']:assert np.array_equal(geom[k],ref[k]),k
 top=pv.read(SURFACE/'exposed_surface.vtp');solid=pv.read(SURFACE/'clot_solid.vtp');vessel=pv.read(SURFACE/'vessel_surface.vtp')
 assert solid.n_open_edges==0 and np.all(top.cell_data['surface_tag']==1)
 tri=top.points[top.faces.reshape(-1,4)[:,1:]];cent=tri.mean(1);cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);areas=np.linalg.norm(cross,axis=1)/2;face_n=cross/(2*areas[:,None]);assert np.all(areas>0)
 owner=ids[cKDTree(X[ids]).query(cent)[1]];mapped=np.bincount(owner,weights=areas,minlength=len(X));area_rel=float(np.max(abs(mapped[ids]-ref['areas'][ids])/ref['areas'][ids]));assert area_rel<2e-6,area_rel
 normals=np.stack([np.bincount(owner,weights=areas*face_n[:,k],minlength=len(X)) for k in range(3)],axis=1)[ids];coherence=np.linalg.norm(normals,axis=1)/mapped[ids];normals/=np.linalg.norm(normals,axis=1)[:,None]
 loc=vtk.vtkStaticCellLocator();loc.SetDataSet(top);loc.BuildLocator();closest=[];nearest=[]
 for x in X[ids]:
  q=[0.,0.,0.];ci=vtk.mutable(0);su=vtk.mutable(0);dist=vtk.mutable(0.);loc.FindClosestPoint(x,q,ci,su,dist);closest.append(q);nearest.append(int(ci))
 closest=np.array(closest);nearest=np.array(nearest);dot=np.einsum('ij,ij->i',normals,face_n[nearest])
 # Closed-solid occupancy independently checks the oriented exposed mesh.
 def inside(mesh,q):return np.asarray(pv.PolyData(q).select_enclosed_points(mesh,tolerance=1e-9,check_surface=True)['SelectedPoints'],bool)
 eps=1e-6;plus=closest+eps*normals;minus=closest-eps*normals;outside=~inside(solid,plus);inward=inside(solid,minus);in_lumen=inside(vessel,plus)
 # Every assigned face is exposed; attachment-neighbour nodes near the rim
 # remain loadable only through their original exposed-face area allocation.
 valid=outside&inward&in_lumen&(dot>0);bad=ids[~valid]+1
 angle=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',normals,ref['tangent'][ids])/np.linalg.norm(ref['tangent'][ids],axis=1),-1,1)))
 np.savez_compressed(O/'REFERENCE_NORMAL_FIELD.npz',atom_ids=ids+1,X_ref=X[ids],normal_ref=normals,tangent_ref=ref['tangent'][ids],area_m2=ref['areas'][ids],mapped_area_m2=mapped[ids],face_owner_atom_ids=owner+1,face_ids=np.arange(len(tri)),face_area_m2=areas,face_normals=face_n,face_centroids_m=cent,closest_surface_point_m=closest,closest_face_ids=nearest,normal_nearest_face_dot=dot,area_normal_coherence=coherence,angle_deg=angle,outside_clot_probe=outside,inside_clot_probe=inward,inside_lumen_probe=in_lumen)
 audit=dict(status='PASS' if not len(bad) else 'FAIL',load_count=len(ids),source_sha256={str(SURFACE/name):sha(SURFACE/name) for name in ['exposed_surface.vtp','clot_solid.vtp','attachment_surface.vtp','vessel_surface.vtp']},geometry_source_matches=True,method='Use oriented triangles of archived exposed surface only; assign each original face once to nearest original loaded material point; normalize area-weighted sum of face normals. Closed original solid defines outward orientation. Attachment faces excluded.',face_count=len(tri),total_original_load_area_m2=float(ref['areas'].sum()),mesh_exposed_area_m2=float(areas.sum()),max_area_relative_difference=area_rel,unit_error=float(abs(np.linalg.norm(normals,axis=1)-1).max()),normal_sign='OUTWARD_FROM_CLOT_INTO_LUMEN',probe_distance_m=eps,outside_plus_count=int(outside.sum()),inside_minus_count=int(inward.sum()),inside_lumen_plus_count=int(in_lumen.sum()),ambiguous_atom_ids=bad,normal_nearest_face_dot_min=float(dot.min()),area_normal_coherence_min=float(coherence.min()),angles_deg_quantiles=np.quantile(angle,[0,.25,.5,.75,1]),edge_rule='Use only owned exposed triangles, no attachment blending; record closest-point sidedness, local-normal alignment and area-vector coherence for every ID; do not replace ambiguous normals.',source_face_indices='0-based indices into archived exposed_surface.vtp',no_geometry_or_area_change=True)
 dump(O/'NORMAL_MAPPING_AUDIT.json',audit);print(json.dumps(audit,default=cv,ensure_ascii=False,indent=2));assert not len(bad),'AMBIGUOUS_NORMALS'
 from external import External,HistoricalExternal
 ns=SimpleNamespace(ref=ref,n=n,N=len(X));old=HistoricalExternal(ns);new=External(ns,query=old.query);checks=[]
 states=[dict(x=X.copy(),velocity=np.zeros_like(X),d=np.zeros(len(ref['bonds'])))]
 key=min((B/'full/keyframes').glob('*.npz'),key=lambda p:abs(float(np.load(p)['time'])-.0005));states.append(load_npz(key))
 for si,s in enumerate(states):
  for t in [0.,.0001,.0005,.0008,.0011605,.00233]:
   u=old.compute(s['x'],s['velocity'],s['d'],t);v=new.compute(s['x'],s['velocity'],s['d'],t);err=float(abs(np.linalg.norm(u['load'],axis=1)-np.linalg.norm(v['load'],axis=1)).max());assert err<1e-20
   assert not v['load'][ref['areas'][:n]==0].any() and np.isfinite(v['load']).all()
   for k in ['drag','numeric','wall_elastic','wall_damping']:assert np.array_equal(u[k],v[k]),k
   assert np.allclose(v['external_force']-u['external_force'],v['load']-u['load'],atol=1e-20,rtol=0)
   checks.append(dict(coordinate_state=si,time_s=t,max_force_magnitude_error_N=err,normal_force_sum_norm_N=float(np.linalg.norm(v['load'],axis=1).sum())))
 dump(O/'LOAD_UNIT_CHECKS.json',dict(status='PASS',checks=checks,nonload_nodes_zero=True,anchor_load_zero=True,original_tangent_norm_range=[float(np.linalg.norm(ref['tangent'][ids],axis=1).min()),float(np.linalg.norm(ref['tangent'][ids],axis=1).max())],original_tangent_magnitude_preserved=True,amp_and_compute_inherited_without_override=True,no_duplicate_tangent_load=True,units='N',nonload_force_channels_unchanged=True))
 (O/'P0_CLEANUP_AUDIT.md').write_text('# P0 code audit\n\nHistorical Native, Integrator, Mechanics, Journal, Checks and TypedLPS are reused by import without edits; External only substitutes the reference direction operand. Full C++ library is reused by SHA; no rebuild. Previous graph_ops and candidate definitions are retained. No old file can be proven both unreferenced and without historical reproduction value, so no historical files were deleted. New task paths are isolated; no git reset/clean.\n')
 dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',dict(status='PASS',baseline_manifest_verified_files=len(manifest['files']),previous_round_verified_files=len(prev['files']),library=str(LIB),library_sha256=sha(LIB),historical_driver_sha256=freeze['driver_sha256'],input_sha256=freeze['input_sha256'],graph_ops_sha256=sha(PREV_A/'graph_ops.py'),Cpp_changed=False,compiled=False))
 dump(O/'PROTECTED_INPUTS.json',{str(p):sha(p) for p in [LIB,B/'PAYLOAD_MANIFEST.json',PREV/'PAYLOAD_MANIFEST.json',R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',Path('/workspace/y_clot/env.sh'),*[HIST_A/name for name in freeze['driver_sha256']],PREV_A/'graph_ops.py']})
 dump(O/'PREPARED.json',dict(status='PASS',wall_s=time.perf_counter()-start,only_physical_change='REFERENCE_OUTWARD_CLOT_SURFACE_NORMAL',normal_sha256=sha(O/'REFERENCE_NORMAL_FIELD.npz')))
if __name__=='__main__':main()
