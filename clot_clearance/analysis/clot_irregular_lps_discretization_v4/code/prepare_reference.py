"""Conservative mapped control cells and deterministic, phase-volume constrained jitter.

Auxiliary cells are quadrature only. LPS connections come from a fresh KD-tree.
No material phase changes, remeshing in current space, or outcome-based selection.
"""
from common import *
from itertools import product,permutations
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import csv,time,shutil
COUNTS=(25,25,13);CORNERS=np.array(list(product([0,1],repeat=3)),int)

class MappedCells:
 def __init__(self):
  shape=tuple(n+1 for n in COUNTS);self.corner_grid=np.indices(shape).reshape(3,-1).T
  self.reference=(self.corner_grid-np.array(COUNTS)/2).astype(float)
  cells=np.indices(COUNTS).reshape(3,-1).T
  self.corner_ids=np.ravel_multi_index((cells[:,None,:]+CORNERS).transpose(2,0,1),shape)
  inside=((self.corner_grid>0)&(self.corner_grid<COUNTS)).all(1);self.inside=np.flatnonzero(inside)
  self.support=np.ravel_multi_index((self.corner_grid[inside,None,:]-1+CORNERS).transpose(2,0,1),COUNTS)
  index={tuple(q):i for i,q in enumerate(CORNERS)};tets=[]
  for perm in permutations(range(3)):
   q=np.zeros(3,int);ids=[index[tuple(q)]]
   for axis in perm:q=q.copy();q[axis]+=1;ids.append(index[tuple(q)])
   x=CORNERS[ids];sign=np.linalg.det((x[1:]-x[0]).T)
   if sign<0:ids[1],ids[2]=ids[2],ids[1]
   tets.append(ids)
  self.pattern=np.array(tets);self.tets=self.corner_ids[:,self.pattern].reshape(-1,4)
 def map(self,u):
  x=self.reference.copy();x[self.inside]+=u[self.support].mean(1);return x
 def volumes(self,x,phase=None):
  q=x[self.tets];a,b,c=(q[:,k]-q[:,0] for k in [1,2,3]);tv=np.einsum('ij,ij->i',a,np.cross(b,c))/6;vol=tv.reshape(-1,6).sum(1)
  if phase is None:return vol,tv
  weights=np.repeat(phase,6);g1=np.cross(b,c)/6;g2=np.cross(c,a)/6;g3=np.cross(a,b)/6;g0=-g1-g2-g3;gc=np.zeros_like(x)
  for j,g in enumerate([g0,g1,g2,g3]):
   for k in range(3):gc[:,k]+=np.bincount(self.tets[:,j],weights=weights*g[:,k],minlength=len(x))
  gu=np.zeros((np.prod(COUNTS),3))
  for j in range(8):
   for k in range(3):gu[:,k]+=np.bincount(self.support[:,j],weights=gc[self.inside,k]/8,minlength=len(gu))
  return vol,tv,gu

def main():
 start=time.perf_counter();(O/'provenance').mkdir(parents=True,exist_ok=True)
 assert (O/'PERFORMANCE_BASELINE_V3.json').exists() and (O/'PERFORMANCE_OPTIMIZATION.json').exists(),'P0/P1 must precede reference changes'
 assert not (O/'provenance/reference.npz').exists(),'Do not replace a frozen reference'
 old=load_npz(V3_O/'provenance/reference.npz');field=load_npz(V3_O/'BIPHASIC_FIELD.npz');gates=json.loads((A/'admission_gates.json').read_text());dump(O/'ADMISSION_GATES.json',gates)
 X0=old['X'];n=int(old['clot_count']);N=len(X0);assert n==8125 and N==10000
 phase=field['material_phase'];strong=phase[:n]==1
 outer=np.any(np.isclose(X0[:n],X0[:n].min(0),rtol=0,atol=1e-14)|np.isclose(X0[:n],X0[:n].max(0),rtol=0,atol=1e-14),axis=1)
 notch=np.isclose(X0[:n,2],0,rtol=0,atol=1e-14)|(np.isclose(X0[:n,2],DX,rtol=0,atol=1e-14)&(X0[:n,0]<=-DX+1e-14))
 fixed=outer|notch
 rng=np.random.default_rng(20261019);direction=rng.normal(size=(n,3));direction/=np.linalg.norm(direction,axis=1)[:,None]
 proposal=.25*rng.random(n)[:,None]**(1/3)*direction;proposal[fixed]=0;u=proposal.copy();cells=MappedCells();trace=[]
 target=float(strong.sum())
 for iteration in range(30):
  corners=cells.map(u);vol,tv,grad=cells.volumes(corners,strong);res=float(vol[strong].sum()-target);trace.append(dict(iteration=iteration,strong_volume_error_dx3=res))
  if abs(res)<1e-10:break
  grad[fixed]=0;step=res/np.sum(grad*grad);u-=step*grad
  length=np.linalg.norm(u,axis=1);u*=np.minimum(1.,.25/np.maximum(length,1e-30))[:,None];u[fixed]=0
 else:raise AssertionError('Phase-volume geometric projection did not converge')
 assert tv.min()/(1/6)>.1 and np.max(np.linalg.norm(u,axis=1))<=.25*(1+1e-13)
 assert abs(vol.sum()/n-1)<1e-12 and abs(vol[strong].sum()/target-1)<1e-12
 X=X0.copy();X[:n]+=u*DX;volumes=np.r_[vol*V,np.full(N-n,V)]
 assert np.array_equal(X[n:],X0[n:]) and np.array_equal(X[:n][fixed],X0[:n][fixed])
 tree=cKDTree(X);minimum=float(tree.query(X,k=2)[0][:,1].min());assert minimum>=.5*DX
 b=tree.query_pairs(HORIZON,output_type='ndarray');b.sort(1);b=b[np.lexsort((b[:,1],b[:,0]))];i,j=b.T;r0=np.linalg.norm(X[j]-X[i],axis=1)
 h=np.minimum(1.,(HORIZON+DX/2-r0)/DX);assert h.min()>0
 m=np.bincount(i,weights=h*r0*volumes[j],minlength=N)+np.bincount(j,weights=h*r0*volumes[i],minlength=N)
 bt=np.where((b<n).all(1),1,np.where((b>=n).all(1),3,2));side=X[:,2]>1e-14;ids=np.flatnonzero((side[i]!=side[j])&(bt==1));alpha=-X[i[ids],2]/(X[j[ids],2]-X[i[ids],2]);hit=X[i[ids]]+alpha[:,None]*(X[j[ids]]-X[i[ids]])
 cut=np.zeros(len(b),bool);cut[ids]=hit[:,0]<=-50e-6+1e-14
 live=b[~cut];nc,labels=connected_components(csr_matrix((np.ones(2*len(live)),(np.r_[live[:,0],live[:,1]],np.r_[live[:,1],live[:,0]])),shape=(N,N)),directed=False);assert nc==1
 near=tree.query_pairs(DX*(1-1e-12),output_type='ndarray');codes=b[:,0]*N+b[:,1];livecodes=codes[~cut];nonbond_contact=near[~np.isin(near[:,0]*N+near[:,1],livecodes)];assert len(nonbond_contact)==0,'Unloaded pre-notch contact introduced'
 ref=dict(X=X,regular_X=X0.copy(),bonds=b,bond_types=bt,original_bond_ids=np.arange(1,len(b)+1),r0=r0,taper=h,m=m,V=volumes,clot_count=n,free_dof=old['free_dof'],initial_cut=cut,atom_ids=old['atom_ids'],control_corner_ids=cells.corner_ids,control_corners_m=corners*DX,control_tet_pattern=cells.pattern,irregular_fixed=fixed)
 np.savez_compressed(O/'provenance/reference.npz',**ref)
 with (O/'provenance/reference.native').open('w') as f:
  f.write(f'LPSVAR_REFERENCE_V2\n{N} {len(b)} {V:.17g}\n')
  for a in range(N):f.write(f'{a+1} {m[a]:.17g} {volumes[a]:.17g}\n')
  for k,((a,c),r,taper,typ) in enumerate(zip(b,r0,h,bt)):f.write(f'{k+1} {a+1} {c+1} {typ} {r:.17g} {taper:.17g} 0 0 0\n')
 with (O/'provenance/body.data').open('w') as f:
  f.write(f'V4 irregular reference, persistent original atom IDs\n\n{N} atoms\n{len(b)} bonds\n\n2 atom types\n3 bond types\n\n')
  for axis in 'xyz':f.write(f'-0.003 0.003 {axis}lo {axis}hi\n')
  f.write(f'\nMasses\n\n1 {RHO*V:.17g}\n2 {RHO*V:.17g}\n\nAtoms # bond\n\n')
  for k,x in enumerate(X):f.write(f'{k+1} 1 {1 if k<n else 2} '+' '.join(f'{v:.17g}' for v in x)+'\n')
  f.write('\nBonds\n\n')
  for k,((a,c),t) in enumerate(zip(b,bt)):f.write(f'{k+1} {t} {a+1} {c+1}\n')
 bulk=bt==1;on=S0[bt].copy();off=SC[bt].copy();on[bulk]=np.where(phase[b[bulk]]==1,.25,.10).mean(1);off[bulk]=np.where(phase[b[bulk]]==1,.30,.12).mean(1)
 props={k:field[k].copy() for k in ['atom_ids','material_phase','strong_domain_id','grid_domain_id','seed_centers_lattice']}
 props.update(bond_ids=np.arange(1,len(b)+1),bond_s0=on,bond_sc=off,cross_phase=bulk&(phase[i]!=phase[j]),bulk_reference_weights=h*r0*volumes[i]*volumes[j]*(1/m[i]+1/m[j]))
 np.savez_compressed(O/'BIPHASIC_FIELD.npz',**props)
 for name in ['BIPHASIC_MATERIAL_SPEC.json','FROZEN_SHEAR_DRIVE.json','SHEAR_DRIVE.npz','MULTIFRAGMENT_OBJECTIVE.json']:shutil.copy2(V3_O/name,O/name)
 np.savez_compressed(O/'REGULAR_TO_IRREGULAR_MAP.npz',atom_ids=old['atom_ids'],regular_X_m=X0,irregular_X_m=X,displacement_m=X-X0,uniform_ball_proposal_m=proposal*DX,phase_volume_projection_correction_m=(u-proposal)*DX)
 oldcodes=old['bonds'][:,0]*N+old['bonds'][:,1];allcodes=np.union1d(oldcodes,codes);oldid=np.searchsorted(oldcodes,allcodes);newid=np.searchsorted(codes,allcodes);oldok=(oldid<len(oldcodes))&(oldcodes[np.minimum(oldid,len(oldcodes)-1)]==allcodes);newok=(newid<len(codes))&(codes[np.minimum(newid,len(codes)-1)]==allcodes)
 with (O/'REFERENCE_BOND_TOPOLOGY.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['atom_i','atom_j','v3_bond_id','v4_bond_id','status'])
  w.writerows((int(c//N+1),int(c%N+1),int(a+1) if ao else -1,int(bid+1) if bo else -1,'RETAINED' if ao and bo else 'ADDED' if bo else 'REMOVED') for c,a,bid,ao,bo in zip(allcodes,oldid,newid,oldok,newok))
 with (O/'NODAL_VOLUME_AUDIT.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['original_atom_id','anchor','phase','volume_m3','mass_kg','volume_over_dx3','displacement_m'])
  w.writerows((k+1,int(k>=n),int(phase[k]),volumes[k],RHO*volumes[k],volumes[k]/V,float(np.linalg.norm(X[k]-X0[k]))) for k in range(N))
 neighbor=np.bincount(b.ravel(),minlength=N);direction=(X[j]-X[i])/r0[:,None];tensor=direction[bulk].T@direction[bulk]/bulk.sum()
 audit=dict(status='REFERENCE_GEOMETRY_PASS',seed=20261019,material_nodes=n,anchor_nodes=N-n,distribution='Uniform isotropic proposal in a 3D ball of radius 12.5 um; fixed exterior/notch nodes; deterministic geometric strong-volume projection with radial hard cap. Accepted displacement distribution is conditioned by these constraints.',displaced_nodes=int((np.linalg.norm(u,axis=1)>1e-13).sum()),fixed_outer_or_notch_nodes=int(fixed.sum()),max_displacement_m=float(np.linalg.norm(u,axis=1).max()*DX),mean_displacement_m=float(np.linalg.norm(u,axis=1).mean()*DX),projection_max_correction_m=float(np.linalg.norm(u-proposal,axis=1).max()*DX),phase_volume_projection=trace,minimum_point_spacing_m=minimum,minimum_tet_signed_volume_ratio=float(tv.min()*6),material_volume_m3=float(volumes[:n].sum()),material_mass_kg=float(RHO*volumes[:n].sum()),center_envelope_volume_m3=float(np.prod(np.ptp(X[:n],axis=0))),strong_volume_fraction=float(volumes[:n][strong].sum()/volumes[:n].sum()),strong_nodes=int(strong.sum()),strong_regions=10,phase_ID_membership_exact=True,nodal_volume_min_max_m3=[float(volumes[:n].min()),float(volumes[:n].max())],original_bonds=len(oldcodes),new_reference_bonds=len(codes),added_bonds=int((newok&~oldok).sum()),removed_bonds=int((oldok&~newok).sum()),initial_cut_bonds=int(cut.sum()),initial_components=nc,initial_unbonded_contact_pairs=len(nonbond_contact),neighbor_count_min_mean_max=[int(neighbor[:n].min()),float(neighbor[:n].mean()),int(neighbor[:n].max())],bulk_direction_tensor=tensor,bulk_direction_tensor_eigenvalues=np.linalg.eigvalsh(tensor),thresholds_unchanged=True,wall_s=time.perf_counter()-start,source_V3_reference_sha256=sha(V3_O/'provenance/reference.npz'),reference_sha256=sha(O/'provenance/reference.npz'),quadrature='Shared mapped control-cell corners; six positive oriented auxiliary tetra per cell, no rescaling of volumes. Total and phase volumes constrained geometrically. Anchors retain separate original volume.',notch='Same z=0 / x<=-50um physical segment cut; z=0 center layer and opposing nearest notch-face nodes fixed to avoid artificial initial contact.')
 dump(O/'IRREGULAR_REFERENCE_GEOMETRY.json',audit);print('IRREGULAR_REFERENCE_PASS',json.dumps(audit,default=cv),flush=True)
if __name__=='__main__':main()
