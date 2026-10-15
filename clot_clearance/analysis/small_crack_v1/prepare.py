from common import *
from scipy.spatial import cKDTree
import shutil,difflib
assert not O.exists(),'Never overwrite an archived case'
for p in [O,F,O/'provenance',O/'source_identity',O/'build',O/'preflight',O/'full']:p.mkdir(parents=True,exist_ok=True)
manifest=json.loads((PREV/'PAYLOAD_MANIFEST.json').read_text())
for n,m in manifest['files'].items():assert sha(R/n)==m['sha256'],n
assert sha(OLDLIB)=='08ccf24110e2cf2d07dbc8c0f31c0f368368b19b5427e6e5de9cc108a5ae115d'
shutil.copy2(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
for p in (PREV/'source_identity').iterdir():
 if p.is_file():shutil.copy2(p,O/'source_identity'/p.name)
shutil.copy2(NATIVE_PREV/'source/src/fix_move.h',O/'source_identity/fix_move.h')
shutil.copy2(NATIVE_PREV/'source/src/BPM/pair_bpm_peri.cpp',O/'source_identity/pair_bpm_peri.cpp')
shutil.copy2(NATIVE_PREV/'source/src/BPM/bond_bpm_peri.cpp',O/'source_identity/legacy_bond_bpm_peri.cpp')
protected=[OLDLIB,R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so',Path('/workspace/y_clot/env.sh'),R/'analysis/softening_feasibility_v1/model.py',PREV/'PAYLOAD_MANIFEST.json']
dump(O/'PROTECTED_BASELINES.json',{str(p):sha(p) for p in protected})
X=np.stack(np.meshgrid(*[(np.arange(n)-(n-1)/2)*DX for n in COUNTS],indexing='ij'),axis=-1).reshape(-1,3)
b=cKDTree(X).query_pairs(HORIZON,output_type='ndarray');b.sort(axis=1);b=b[np.lexsort((b[:,1],b[:,0]))];r0=np.linalg.norm(X[b[:,1]]-X[b[:,0]],axis=1)
taper=np.where(abs(r0-HORIZON)<=DX/2,-r0/DX+1+(HORIZON-DX/2)/DX,1.)
m=np.bincount(b.ravel(),weights=np.repeat(taper*r0*V,2),minlength=len(X))
cross=X[b[:,0],2]*X[b[:,1],2]<0;ci=np.flatnonzero(cross);frac=-X[b[ci,0],2]/(X[b[ci,1],2]-X[b[ci,0],2]);hit=X[b[ci,0]]+frac[:,None]*(X[b[ci,1]]-X[b[ci,0]])
cut=np.zeros(len(b),bool);cut[ci]=hit[:,0]<=TIP+1e-14
bottom=np.flatnonzero(X[:,2]==X[:,2].min())+1;top=np.flatnonzero(X[:,2]==X[:,2].max())+1;grip=np.r_[bottom,top];free=np.setdiff1d(np.arange(1,len(X)+1),grip)
free_dof=np.ones(X.shape,bool);free_dof[grip-1,2]=False
node_regions=classify(X);bond_regions=classify(.5*(X[b[:,0]]+X[b[:,1]]))
# Exact segment-plane intersections, including every long-range reference bond.
ix=np.clip(np.floor((hit[:,0]+SIZE[0]/2)/DX).astype(int),0,COUNTS[0]-1);iy=np.clip(np.floor((hit[:,1]+SIZE[1]/2)/DX).astype(int),0,COUNTS[1]-1);cells=ix*COUNTS[1]+iy
face_lower=[];face_upper=[]
for x in np.unique(X[:,0]):
 for y in np.unique(X[:,1]):
  same=np.flatnonzero((X[:,0]==x)&(X[:,1]==y));face_lower.append(same[np.argmax(np.where(X[same,2]<0,X[same,2],-np.inf))]);face_upper.append(same[np.argmin(np.where(X[same,2]>0,X[same,2],np.inf))])
ref=dict(X=X,bonds=b,r0=r0,taper=taper,m=m,V=np.full(len(X),V),atom_ids=np.arange(1,len(X)+1),bond_ids=np.arange(1,len(b)+1),initial_cut=cut,cross=cross,cross_ids=ci,cross_intersections=hit,cross_cells=cells,node_regions=node_regions,bond_regions=bond_regions,face_lower=np.array(face_lower),face_upper=np.array(face_upper),free_dof=free_dof)
np.savez_compressed(O/'provenance/reference.npz',**ref);np.savez_compressed(O/'provenance/boundary_ids.npz',bottom_ids=bottom,top_ids=top,free_ids=free,free_dof=free_dof);np.save(O/'provenance/original_bond_ids.npy',np.arange(1,len(b)+1))
out=['LPSVAR_REFERENCE_V1',f'{len(X)} {len(b)} {V:.17g}']+[f'{i+1} {z:.17g}' for i,z in enumerate(m)]+[f'{k+1} {i+1} {j+1} 1 {r:.17g} {h:.17g} 0 0 0' for k,((i,j),r,h) in enumerate(zip(b,r0,taper))]
(O/'provenance/reference.native').write_text('\n'.join(out)+'\n')
out=['3D single-edge precrack; full immutable reference graph before initial cut','',f'{len(X)} atoms',f'{len(b)} bonds','','1 atom types','1 bond types','']+[f'{-q:.17g} {q:.17g} {axis}lo {axis}hi' for q,axis in zip(BOX,'xyz')]+['','Masses','',f'1 {RHO*V:.17g}','','Atoms # bond','']+[f'{i+1} 1 1 '+' '.join(f'{z:.17g}' for z in x) for i,x in enumerate(X)]+['','Bonds','']+[f'{k+1} 1 {i+1} {j+1}' for k,(i,j) in enumerate(b)]
(O/'provenance/body.data').write_text('\n'.join(out)+'\n')
audit=dict(status='GEOMETRY_GENERATED_PENDING_NATIVE_CHECK',shape='14 x 10 x 22 cell-centered cubic lattice',nodes=len(X),original_reference_bonds=len(b),initial_cut_bonds=int(cut.sum()),initial_live_bonds=int((~cut).sum()),initial_cross_plane_bonds=len(ci),initial_live_cross_plane_bonds=int(np.sum(cross&~cut)),precut_original_bond_ids=np.flatnonzero(cut)+1,initial_live_cross_bond_ids=np.flatnonzero(cross&~cut)+1,physical_voxel_dimensions_m=SIZE,node_coordinate_spans_m=np.ptp(X,axis=0),node_xyz_min_m=X.min(0),node_xyz_max_m=X.max(0),notch_plane='z=0; cut if the ORIGINAL segment-plane intersection has x<=-50 micrometers',front_x_m=TIP,front_y_interval_m=[-SIZE[1]/2,SIZE[1]/2],initial_projected_cut_area_m2=(TIP+SIZE[0]/2)*SIZE[1],remaining_ligament_m=SIZE[0]/2-TIP,full_reference_m_policy='m_i=sum_full_original_neighbors taper*r0*V including initial cuts; never renormalize after precut or natural damage',initial_state_policy='Import kappa=sc,d=1 ONLY on initial cuts before attaching dynamic coordinator; Diss=0, natural event count=0. Keep original records; native topology removes initial cuts.',mass_per_node_kg=RHO*V,normal_control_nodes=len(grip),fully_free_nodes=len(free),controlled_DOFs=len(grip),free_DOFs=int(free_dof.sum()),normal_control_axis='z only; x,y free for ALL nodes',horizon_scale=dict(width=SIZE[0]/HORIZON,thickness=SIZE[1]/HORIZON,height=SIZE[2]/HORIZON,ligament=(SIZE[0]/2-TIP)/HORIZON,tip_to_control_plane=(SIZE[2]-DX)/2/HORIZON,thickness_core_m=SIZE[1]-2*HORIZON),region_labels=REGIONS,region_rule='reference midpoint/node center; notch within radial delta first (split side/core), crack flank x<tip and |z|<=delta/2, actual outermost grip-plane distance<=delta, other external surface distance<=delta, remaining interior',physical_limits='Finite horizon and free thickness boundaries remain. Not plane strain, not a continuum-size-converged fracture experiment.',no_parameter_or_geometry_scan=True)
dump(O/'GEOMETRY_AND_INITIAL_CUT_AUDIT.json',audit)
print('GEOMETRY',len(X),len(b),'CUT',cut.sum(),'LIVE_CROSS',np.sum(cross&~cut))
