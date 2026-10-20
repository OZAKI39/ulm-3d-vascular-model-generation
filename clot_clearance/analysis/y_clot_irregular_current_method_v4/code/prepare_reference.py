"""V4 quadrature on the original Y material mask; no new shape or notch."""
from common import *
from itertools import product,permutations
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import csv,time,shutil

BITS=np.array(list(product([0,1],repeat=3)),int)
class MappedCells:
 def __init__(self,grid):
  self.n=len(grid);lookup={tuple(q):i for i,q in enumerate(grid)}
  keys=(2*grid[:,None,:]+2*BITS-1).reshape(-1,3)
  unique,inv=np.unique(keys,axis=0,return_inverse=True);self.reference=unique/2.;self.corner_ids=inv.reshape(-1,8)
  support=np.array([[lookup.get(tuple((q+2*b-1)//2),-1) for b in BITS] for q in unique])
  self.inside=np.flatnonzero((support>=0).all(1));self.support=support[self.inside]
  index={tuple(q):i for i,q in enumerate(BITS)};pattern=[]
  for perm in permutations(range(3)):
   q=np.zeros(3,int);ids=[index[tuple(q)]]
   for axis in perm:q=q.copy();q[axis]+=1;ids.append(index[tuple(q)])
   if np.linalg.det((BITS[ids][1:]-BITS[ids][0]).T)<0:ids[1],ids[2]=ids[2],ids[1]
   pattern.append(ids)
  self.pattern=np.array(pattern);self.tets=self.corner_ids[:,self.pattern].reshape(-1,4)
 def map(self,u):
  corners=self.reference.copy();corners[self.inside]+=u[self.support].mean(1);return corners
 def volumes(self,corners,phase):
  q=corners[self.tets];a,b,c=(q[:,k]-q[:,0] for k in [1,2,3]);tv=np.einsum('ij,ij->i',a,np.cross(b,c))/6;vol=tv.reshape(-1,6).sum(1)
  weights=np.repeat(phase,6);g1=np.cross(b,c)/6;g2=np.cross(c,a)/6;g3=np.cross(a,b)/6;g0=-g1-g2-g3;gc=np.zeros_like(corners)
  for j,g in enumerate([g0,g1,g2,g3]):
   for k in range(3):gc[:,k]+=np.bincount(self.tets[:,j],weights=weights*g[:,k],minlength=len(corners))
  gu=np.zeros((self.n,3))
  for j in range(8):
   for k in range(3):gu[:,k]+=np.bincount(self.support[:,j],weights=gc[self.inside,k]/8,minlength=self.n)
  return vol,tv,gu

def main():
 start=time.perf_counter();assert not (O/'provenance/reference.npz').exists()
 for path in [O/'provenance',O/'build',O/'checks',F/'figures',F/'videos']:path.mkdir(parents=True,exist_ok=True)
 old=load_npz(Y_OLD/'provenance/reference.npz');field=load_npz(Y_OLD/'BIPHASIC_FIELD.npz')
 X0=old['X'];n=int(old['clot_count']);N=len(X0);assert n==1927 and N==3488 and not old['initial_cut'].any()
 grid=np.rint(X0[:n]/DX).astype(int);assert abs(grid*DX-X0[:n]).max()<1e-14
 lookup={tuple(q):i for i,q in enumerate(grid)}
 offsets=np.array(list(product([-1,0,1],repeat=3)));fixed=np.array([any(tuple(q+e) not in lookup for e in offsets) for q in grid])
 # Original wall/interface material nodes stay at the historical attachment.
 fixed[np.unique(old['bonds'][old['bond_types']==2,0])]=True
 phase=field['material_phase'];strong=phase[:n]==1;assert strong.sum()==882
 rng=np.random.default_rng(20261019);direction=rng.normal(size=(n,3));direction/=np.linalg.norm(direction,axis=1)[:,None]
 proposal=.25*rng.random(n)[:,None]**(1/3)*direction;proposal[fixed]=0;u=proposal.copy();cells=MappedCells(grid);trace=[]
 for iteration in range(30):
  corners=cells.map(u);vol,tv,grad=cells.volumes(corners,strong);res=float(vol[strong].sum()-strong.sum());trace.append(dict(iteration=iteration,strong_volume_error_dx3=res))
  if abs(res)<1e-10:break
  grad[fixed]=0;assert np.sum(grad*grad)>0;u-=res/np.sum(grad*grad)*grad
  length=np.linalg.norm(u,axis=1);u*=np.minimum(1.,.25/np.maximum(length,1e-30))[:,None];u[fixed]=0
 else:raise RuntimeError('Phase volume projection failed before any dynamics')
 assert tv.min()*6>.1 and abs(vol.sum()/n-1)<1e-12 and abs(vol[strong].sum()/strong.sum()-1)<1e-12
 X=X0.copy();X[:n]+=u*DX;volumes=np.r_[vol*V,old['V'][n:]]
 assert np.linalg.norm(X-X0,axis=1).max()<=12.5e-6*(1+1e-12)
 from pd_clot.lammps_fragmentation.native_y_wall_query import NativeYWallQuery
 from pd_clot.lammps_fragmentation.y_mural_clot_geometry import locator,hit,make_surface,inside
 from pd_clot.lammps_fragmentation.y_clot_morphology_v3 import profile
 wall=NativeYWallQuery();configpath=R/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/V3_GEOMETRY_CONFIG.json'
 config=json.loads(configpath.read_text());solid,top,attachment=make_surface(wall.sources,config,profile_fn=profile)
 assert inside(solid,X[:n]).all(),'Material center moved outside historical continuous clot'
 sl=locator(solid);tl=locator(top);al=locator(attachment)
 tree=cKDTree(X);minimum=float(tree.query(X,k=2)[0][:,1].min());assert minimum>=.5*DX
 candidates=tree.query_pairs(HORIZON*(1+1e-12),output_type='ndarray');candidates.sort(1)
 oldcodes=old['bonds'][:,0]*N+old['bonds'][:,1];oldset=set(oldcodes.tolist());moved=np.linalg.norm(X-X0,axis=1)>1e-15
 keep=[];rejected={1:0,2:0,3:0};ray_queries=0
 for i,j in candidates:
  kind=1 if j<n else 2 if i<n else 3
  if not (moved[i] or moved[j]):valid=int(i*N+j) in oldset
  elif kind==1:valid=not hit(sl,X[i],X[j])[0];ray_queries+=1
  elif kind==2:valid=hit(al,X[i],X[j])[0] and not hit(tl,X[i],X[j])[0];ray_queries+=1
  else:raise AssertionError('Fixed anchors must not move')
  if valid:keep.append((i,j))
  else:rejected[kind]+=1
 b=np.asarray(keep,int);b=b[np.lexsort((b[:,1],b[:,0]))];i,j=b.T
 r0=np.linalg.norm(X[j]-X[i],axis=1);h=np.minimum(1.,(HORIZON+DX/2-r0)/DX)
 m=np.bincount(i,weights=h*r0*volumes[j],minlength=N)+np.bincount(j,weights=h*r0*volumes[i],minlength=N)
 bt=np.where((b<n).all(1),1,np.where((b>=n).all(1),3,2));cut=np.zeros(len(b),bool);codes=i*N+j
 assert np.array_equal(b[bt==3],old['bonds'][old['bond_types']==3])
 nc,labels=connected_components(csr_matrix((np.ones(2*len(b)),(np.r_[i,j],np.r_[j,i])),shape=(N,N)),directed=False);assert nc==1
 near=tree.query_pairs(DX*(1-1e-12),output_type='ndarray');nonbond=near[~np.isin(near[:,0]*N+near[:,1],codes)];assert len(nonbond)==0,'Initial new unbonded contact introduced'
 ref=dict(X=X,regular_X=X0.copy(),bonds=b,bond_types=bt,original_bond_ids=np.arange(1,len(b)+1),r0=r0,taper=h,m=m,V=volumes,clot_count=n,free_dof=old['free_dof'],initial_cut=cut,atom_ids=np.arange(1,N+1),control_corner_ids=cells.corner_ids,control_corners_m=corners*DX,control_tet_pattern=cells.pattern,irregular_fixed=fixed)
 np.savez_compressed(O/'provenance/reference.npz',**ref)
 with (O/'provenance/reference.native').open('w') as f:
  f.write(f'LPSVAR_REFERENCE_V2\n{N} {len(b)} {V:.17g}\n')
  for a in range(N):f.write(f'{a+1} {m[a]:.17g} {volumes[a]:.17g}\n')
  for k,((a,c),r,t,typ) in enumerate(zip(b,r0,h,bt)):f.write(f'{k+1} {a+1} {c+1} {typ} {r:.17g} {t:.17g} 0 0 0\n')
 with (O/'provenance/body.data').open('w') as f:
  f.write(f'Original Y geometry with V4 irregular reference\n\n{N} atoms\n{len(b)} bonds\n\n2 atom types\n3 bond types\n\n')
  for axis in 'xyz':f.write(f'-0.03 0.03 {axis}lo {axis}hi\n')
  f.write(f'\nMasses\n\n1 {RHO*V:.17g}\n2 {RHO*V:.17g}\n\nAtoms # bond\n\n')
  for k,x in enumerate(X):f.write(f'{k+1} 1 {1 if k<n else 2} '+' '.join(f'{v:.17g}' for v in x)+'\n')
  f.write('\nBonds\n\n')
  for k,((a,c),t) in enumerate(zip(b,bt)):f.write(f'{k+1} {t} {a+1} {c+1}\n')
 bulk=bt==1;on=S0[bt].copy();off=SC[bt].copy();on[bulk]=np.where(phase[b[bulk]]==1,.25,.10).mean(1);off[bulk]=np.where(phase[b[bulk]]==1,.30,.12).mean(1)
 props={k:v.copy() for k,v in field.items() if k not in ['bond_ids','bond_s0','bond_sc']}
 props.update(bond_ids=np.arange(1,len(b)+1),bond_s0=on,bond_sc=off,cross_phase=bulk&(phase[i]!=phase[j]));np.savez_compressed(O/'BIPHASIC_FIELD.npz',**props)
 for name in ['BIPHASIC_MATERIAL_SPEC.json','FROZEN_SHEAR_DRIVE.json']:shutil.copy2(Y_OLD/name,O/name)
 shutil.copy2(A/'admission_gates.json',O/'ADMISSION_GATES.json')
 np.savez_compressed(O/'REGULAR_TO_IRREGULAR_MAP.npz',atom_ids=np.arange(1,N+1),regular_X_m=X0,irregular_X_m=X,displacement_m=X-X0,uniform_ball_proposal_m=proposal*DX,phase_volume_projection_correction_m=(u-proposal)*DX)
 oldmap={int(c):k+1 for k,c in enumerate(oldcodes)};newmap={int(c):k+1 for k,c in enumerate(codes)}
 with (O/'REFERENCE_BOND_TOPOLOGY.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['atom_i','atom_j','old_y_bond_id','v4_y_bond_id','status'])
  for c in np.union1d(oldcodes,codes):w.writerow([c//N+1,c%N+1,oldmap.get(int(c),-1),newmap.get(int(c),-1),'RETAINED' if c in oldmap and c in newmap else 'ADDED' if c in newmap else 'REMOVED'])
 with (O/'NODAL_VOLUME_AUDIT.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['original_atom_id','anchor','phase','volume_m3','mass_kg','volume_over_dx3','displacement_m'])
  for k in range(N):w.writerow([k+1,int(k>=n),int(phase[k]),volumes[k],RHO*volumes[k],volumes[k]/V,float(np.linalg.norm(X[k]-X0[k]))])
 neighbor=np.bincount(b.ravel(),minlength=N)
 audit=dict(status='REFERENCE_GEOMETRY_PASS',material_nodes=n,anchor_nodes=N-n,displaced_nodes=int(moved[:n].sum()),fixed_boundary_and_interface_nodes=int(fixed.sum()),seed=20261019,max_displacement_m=float(np.linalg.norm(X-X0,axis=1).max()),minimum_point_spacing_m=minimum,minimum_control_tet_volume_ratio=float(tv.min()*6),material_volume_m3=float(volumes[:n].sum()),material_mass_kg=float(RHO*volumes[:n].sum()),strong_nodes=int(strong.sum()),strong_volume_fraction=float(volumes[:n][strong].sum()/volumes[:n].sum()),retained_original_domain_labels=np.unique(field['strong_domain_id'][:n][strong]),nodal_volume_min_max_m3=[float(volumes[:n].min()),float(volumes[:n].max())],original_bonds=len(oldcodes),new_reference_bonds=len(codes),added_bonds=len(set(codes)-oldset),removed_bonds=len(oldset-set(codes)),bond_type_counts={str(k):int((bt==k).sum()) for k in [1,2,3]},initial_components=nc,initial_cut_bonds=0,initial_unbonded_contact_pairs=len(nonbond),neighbor_count_min_mean_max=[int(neighbor[:n].min()),float(neighbor[:n].mean()),int(neighbor[:n].max())],phase_volume_projection=trace,void_filter_rejected=rejected,changed_reference_ray_queries=ray_queries,original_boundary_positions_exact=True,anchor_positions_exact=True,reference_sha256=sha(O/'provenance/reference.npz'),source_Y_reference_sha256=sha(Y_OLD/'provenance/reference.npz'),historical_shape_config_sha256=sha(configpath),quadrature='Shared mapped dual cells over the historical discrete material mask. Exterior corners fixed; eight-material-cell interior corners follow average jitter. Six positive tetra per cell. Conserves historical 0.240875 mm3 integration volume, not the separate continuous shape volume.',distribution='Seeded isotropic uniform-ball proposals, boundary/interface constrained; small geometric projection preserves strong-phase integral without rescaling Vi.',wall_s=time.perf_counter()-start)
 dump(O/'IRREGULAR_REFERENCE_GEOMETRY.json',audit)
 print('Y_V4_REFERENCE_PASS',json.dumps(audit,default=cv),flush=True)
if __name__=='__main__':main()
