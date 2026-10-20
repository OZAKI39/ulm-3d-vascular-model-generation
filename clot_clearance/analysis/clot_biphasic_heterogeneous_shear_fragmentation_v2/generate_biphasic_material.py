"""P0 only: original geometry and deterministic connected binary regions, no dynamics."""
from common import *
from scipy.spatial import cKDTree
from scipy.ndimage import label,generate_binary_structure
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import csv,heapq,itertools,subprocess,datetime,shutil,time

def writecsv(path,rows):
 with path.open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def regions():
 rng=np.random.default_rng(SEED);shape=COUNTS;grid=np.indices(shape).transpose(1,2,3,0);phase=np.zeros(shape,np.int16)
 # Nine lateral columns with one selected vertical cell each, and the other
 # central cell: ten candidates. Weak separator sheets remain ordinary material.
 cells=[np.arange(0,8),np.arange(9,16),np.arange(17,25)];zs=[np.arange(0,6),np.arange(7,13)]
 candidates=[(i,j,int(rng.integers(2))) for i in range(3) for j in range(3)]
 central=next(k for i,j,k in candidates if (i,j)==(1,1));candidates.append((1,1,1-central))
 total=round(np.prod(shape)/3);targets=np.full(10,total//10,int);targets[:total%10]+=1
 centers=[]
 for cid,((i,j,k),target) in enumerate(zip(candidates,targets),1):
  allowed=np.zeros(shape,bool);allowed[np.ix_(cells[i],cells[j],zs[k])]=True
  center=np.array([cells[i].mean(),cells[j].mean(),zs[k].mean()]);seed=tuple(np.rint(center).astype(int));centers.append(center)
  anis=rng.uniform(.85,1.15,3);delta=(grid-center)/np.array([3.6,3.6,2.7]);noise=rng.random(shape)
  priority=np.sum((delta*anis)**2,axis=-1)+.25*np.sin(1.3*grid[...,0]+.7*grid[...,1]+rng.uniform(0,6))+.12*noise
  queue=[(float(priority[seed]),seed)];seen={seed};grown=0
  while grown<target:
   _,q=heapq.heappop(queue);phase[q]=cid;grown+=1
   for axis in range(3):
    for sign in [-1,1]:
     z=list(q);z[axis]+=sign;z=tuple(z)
     if all(0<=z[c]<shape[c] for c in range(3)) and allowed[z] and z not in seen:
      seen.add(z);heapq.heappush(queue,(float(priority[z]),z))
 return phase,np.array(centers),candidates

def main():
 started=time.perf_counter();assert not (O/'GEOMETRY_FREEZE.json').exists()
 for p in [O/'provenance',O/'checks',O/'build',F/'figures',F/'videos']:p.mkdir(parents=True,exist_ok=True)
 # Reuse the existing Cartesian/cKDTree/taper/m construction, enlarged only.
 Xm=np.stack(np.meshgrid(*[(np.arange(n)-(n-1)/2)*DX for n in COUNTS],indexing='ij'),axis=-1).reshape(-1,3);n=len(Xm)
 xa=np.unique(Xm[:,0]);ya=np.unique(Xm[:,1]);za=Xm[:,2].min()-DX*np.arange(1,4)
 anchors=np.stack(np.meshgrid(xa,ya,za,indexing='ij'),axis=-1).reshape(-1,3);X=np.vstack([Xm,anchors]);N=len(X)
 b=cKDTree(X).query_pairs(HORIZON,output_type='ndarray');b.sort(axis=1);b=b[np.lexsort((b[:,1],b[:,0]))];r0=np.linalg.norm(X[b[:,1]]-X[b[:,0]],axis=1)
 taper=np.where(abs(r0-HORIZON)<=DX/2,-r0/DX+1+(HORIZON-DX/2)/DX,1.)
 m=np.bincount(b.ravel(),weights=np.repeat(taper*r0*V,2),minlength=N);bt=np.where((b<n).all(1),1,np.where((b>=n).all(1),3,2))
 # Retain the historical single-edge notch, explicitly assign z=0 nodes to
 # the lower side for this odd-sized lattice; no phase-dependent cutting.
 side=X[:,2]>1e-14;ci=np.flatnonzero((side[b[:,0]]!=side[b[:,1]])&(bt==1));frac=-X[b[ci,0],2]/(X[b[ci,1],2]-X[b[ci,0],2]);hit=X[b[ci,0]]+frac[:,None]*(X[b[ci,1]]-X[b[ci,0]])
 cut=np.zeros(len(b),bool);cut[ci]=hit[:,0]<=-50e-6+1e-14
 free=np.zeros(X.shape,bool);free[:n]=True
 ref=dict(X=X,bonds=b,bond_types=bt,original_bond_ids=np.arange(1,len(b)+1),r0=r0,taper=taper,m=m,V=np.full(N,V),clot_count=n,free_dof=free,initial_cut=cut,atom_ids=np.arange(1,N+1))
 np.savez_compressed(O/'provenance/reference.npz',**ref)
 with (O/'provenance/reference.native').open('w') as f:
  f.write(f'LPSVAR_REFERENCE_V1\n{N} {len(b)} {V:.17g}\n')
  for i,v in enumerate(m):f.write(f'{i+1} {v:.17g}\n')
  for k,((i,j),r,h,t) in enumerate(zip(b,r0,taper,bt)):f.write(f'{k+1} {i+1} {j+1} {t} {r:.17g} {h:.17g} 0 0 0\n')
 with (O/'provenance/body.data').open('w') as f:
  f.write(f'Enlarged common specimen; phase-independent original full graph\n\n{N} atoms\n{len(b)} bonds\n\n2 atom types\n3 bond types\n\n')
  for axis in 'xyz':f.write(f'-0.003 0.003 {axis}lo {axis}hi\n')
  f.write(f'\nMasses\n\n1 {RHO*V:.17g}\n2 {RHO*V:.17g}\n\nAtoms # bond\n\n')
  for i,x in enumerate(X):f.write(f'{i+1} 1 {1 if i<n else 2} '+' '.join(f'{v:.17g}' for v in x)+'\n')
  f.write('\nBonds\n\n')
  for k,((i,j),t) in enumerate(zip(b,bt)):f.write(f'{k+1} {t} {i+1} {j+1}\n')
 phase_grid,centers,candidates=regions();domain=phase_grid.ravel();strong=domain>0;phase=np.r_[strong.astype(np.int8),np.full(N-n,-1,np.int8)]
 assert np.array_equal(phase_grid,regions()[0])
 face=generate_binary_structure(3,1);spatial_labels,ndom=label(phase_grid>0,face);_,nweak=label(phase_grid==0,face)
 assert ndom==10 and nweak==1 and .30<=strong.mean()<=.36
 live=b[~cut];graph=csr_matrix((np.ones(len(live)*2),(np.r_[live[:,0],live[:,1]],np.r_[live[:,1],live[:,0]])),shape=(N,N));ncomp,_=connected_components(graph,directed=False);assert ncomp==1
 rows=[];distances=[]
 for c in range(1,11):
  ids=np.flatnonzero(domain==c);points=Xm[ids];volume=len(ids)*V;others=Xm[strong&(domain!=c)];minimum=float(cKDTree(others).query(points)[0].min());distances.append(minimum)
  assert label((phase_grid==c),face)[1]==1
  rows.append(dict(domain_id=c,nodes=len(ids),reference_volume_m3=volume,equivalent_volume_diameter_m=(6*volume/np.pi)**(1/3),center_x_m=points[:,0].mean(),center_y_m=points[:,1].mean(),center_z_m=points[:,2].mean(),x_span_m=np.ptp(points[:,0]),y_span_m=np.ptp(points[:,1]),z_span_m=np.ptp(points[:,2]),minimum_other_domain_center_distance_m=minimum))
 bulk=bt==1;cross=bulk&(phase[b[:,0]]!=phase[b[:,1]]);node_s0=np.where(phase==1,.25,.10);node_sc=np.where(phase==1,.30,.12);on=S0[bt].copy();off=SC[bt].copy();on[bulk]=node_s0[b[bulk]].mean(1);off[bulk]=node_sc[b[bulk]].mean(1)
 weights=taper*r0*V**2*(1/m[b[:,0]]+1/m[b[:,1]])
 np.savez_compressed(O/'BIPHASIC_FIELD.npz',atom_ids=np.arange(1,N+1),material_phase=phase,strong_domain_id=np.r_[domain,np.zeros(N-n,int)],grid_domain_id=phase_grid,seed_centers_lattice=centers,bond_ids=np.arange(1,len(b)+1),bond_s0=on,bond_sc=off,cross_phase=cross,bulk_reference_weights=weights)
 writecsv(O/'STRONG_DOMAINS.csv',rows)
 writecsv(O/'NODE_MATERIAL_PHASE.csv',[dict(original_atom_id=i+1,x_m=x[0],y_m=x[1],z_m=x[2],material_phase='STRONG' if strong[i] else 'WEAK',strong_domain_id=int(domain[i]),s0=node_s0[i],sc=node_sc[i],reference_volume_m3=V) for i,x in enumerate(Xm)])
 writecsv(O/'PAIR_MATERIAL_PROPERTIES.csv',[dict(original_bond_id=i+1,atom_i=int(p[0]+1),atom_j=int(p[1]+1),type=int(bt[i]),s0=on[i],sc=off[i],cross_phase=bool(cross[i]),initial_precrack=bool(cut[i]),reference_weight=weights[i]) for i,p in enumerate(b)])
 audit=dict(status='PASS',material_nodes=n,anchor_nodes=N-n,original_bulk_bonds=int(bulk.sum()),original_interface_bonds=int((bt==2).sum()),anchor_scaffold_bonds=int((bt==3).sum()),original_bonds=len(b),initial_precrack=int(cut.sum()),initial_live_material_and_anchor_components=ncomp,strong_nodes=int(strong.sum()),weak_nodes=int((~strong).sum()),strong_volume_fraction=float(strong.mean()),strong_domain_count=int(ndom),weak_matrix_components=int(nweak),spatial_component_definition='6-neighbor reference-grid face adjacency, not nonlocal material connectivity',domains=rows,minimum_interdomain_node_distance_m=min(distances),cross_phase_material_bonds=int(cross.sum()),cross_phase_fraction_of_bulk=float(cross.sum()/bulk.sum()),cross_phase_initially_live=int((cross&~cut).sum()),strong_domain_ids_preserved=True,phase_dependent_cuts=0,grid_pitch_m=DX,horizon_m=HORIZON,horizon_over_dx=HORIZON/DX,mean_node_s0=float(node_s0[:n].mean()),mean_node_sc=float(node_sc[:n].mean()),weighted_bulk_s0=float(np.average(on[bulk],weights=weights[bulk])),weighted_bulk_sc=float(np.average(off[bulk],weights=weights[bulk])),structure_scale='Constrained 200–400 um-scale grid regions, irregular connected growth. Actual center spans and volume-equivalent diameters reported separately; these are not target fragment identities.',field_sha256=sha(O/'BIPHASIC_FIELD.npz'),wall_s=time.perf_counter()-started)
 dump(O/'PREFLIGHT_STRUCTURE_AUDIT.json',audit)
 dump(O/'BIPHASIC_MATERIAL_SPEC.json',dict(label=HET_LABEL,seed=SEED,generator='One random choice of vertical cells across nine fixed lateral columns plus opposite central cell; seeded anisotropic face-connected priority growth inside separated cells. No dynamics or outcome selection.',candidate_centers=10,strong_fraction_target=1/3,accepted_fraction=[.30,.36],no_transition_band=True,weak=dict(s0=.10,sc=.12),mixed=dict(s0=.175,sc=.210),strong=dict(s0=.25,sc=.30),uniform=dict(s0=.15,sc=.18),G=G,K=K,rho=RHO,interface=dict(s0=S0[2],sc=SC[2]),rule='Arithmetic endpoint thresholds for bulk, unchanged typed interface/scaffold',field_sha256=sha(O/'BIPHASIC_FIELD.npz'),initial_history='Same kappa=0.30,d=1 on original bulk precut in BOTH U/H, kappa=d=0 otherwise. Common failed sentinel exceeds all bulk completion thresholds; no initial physical strain or dissipation.'))
 dump(O/'GEOMETRY_FREEZE.json',dict(status='FROZEN_BEFORE_NATIVE',counts=COUNTS,material_nodes=n,anchor_layers=3,anchor_nodes=N-n,material_center_spans_m=np.ptp(Xm,axis=0),nodal_quadrature_volume_m3=n*V,nominal_envelope_volume_m3=float(np.prod(np.ptp(Xm,axis=0))),volume_note='Preserve existing V=dx^3 even on inclusive end nodes; quadrature volume differs from center-envelope volume.',initial_notch='Historical single edge: z=0 and segment intersection x<=-50 um. New odd grid assigns z=0 nodes to lower side, crossing is (z>0) != (z>0). No cuts depend on phase or domain.',reference_files={name:sha(O/'provenance'/name) for name in ['reference.npz','reference.native','body.data']},field_sha256=sha(O/'BIPHASIC_FIELD.npz'),same_geometry_for_both_cases=True,reference_m='All original neighbors before any cut, same taper*r0*V rule; never recomputed after damage',boundary='Three full bottom anchor layers fixed; all material nodes integrated, no top grip',box='m m m, initial +/-3mm, no fixed physical confinement'))
 dump(O/'P0_SCOPE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),git_sha=subprocess.check_output(['git','-c',f'safe.directory={R}','rev-parse','HEAD'],cwd=R,text=True).strip(),deleted_files=[],reason='No verified disposable unreferenced file; isolate new round',baseline_library_sha256=sha(OLDLIB),baseline_task=str(PREVIOUS),source_sha256={str(p.relative_to(R)):sha(p) for p in [A/'common.py',A/'generate_biphasic_material.py',PREVIOUS_A/'native.py',PREVIOUS_A/'current_surface.py',R/'analysis/softening_feasibility_v1/model.py']}))
 (O/'provenance/git_existing.diff').write_bytes(subprocess.check_output(['git','-c',f'safe.directory={R}','diff'],cwd=R));shutil.copy2(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',O/'provenance/ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
 print('P0_STRUCTURE_PASS',audit,flush=True)
if __name__=='__main__':main()
