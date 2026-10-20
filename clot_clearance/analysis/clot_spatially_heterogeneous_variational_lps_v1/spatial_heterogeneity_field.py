"""One immutable reference-coordinate realization; no outcome-based seeds."""
from common import *
import gstools as gs
from scipy.optimize import brentq
from scipy.ndimage import maximum_filter,minimum_filter,label
import csv
SEED=20261017;LC=300e-6;ALPHA=.30
SPEC=dict(label=HET_LABEL,seed=SEED,nominal_correlation_length_m=LC,alpha=ALPHA,
    covariance='Gaussian C_z(r)=exp(-(r/LC)^2); GSTools rescale=1 explicitly',variance_z=1.,isotropic=True,mode_no=4096,
    bounded_transform='q=tanh(z-c), c solves reference material-volume weighted mean(q)=0',bounds=[-1,1],
    mapping='q_ij=(q_i+q_j)/2; Delta_s=(0.18-0.15)*exp(alpha*q_ij)/Z; sc=0.15+Delta_s; interface unchanged',
    normalization='Z=sum(w exp(alpha*q_ij))/sum(w) over ALL original bulk bonds including original precut; w=taper*r0*V^2*(1/m_i+1/m_j)',
    meaning='Weights proportional to original intact hydrostatic pair stiffness times r0^2. Mean softening interval preserved, not all-path dissipation or Gc.',
    experimental_calibration=False,original_topology_changed=False)
def generate(ref,order=None):
 X=ref['X'];n=int(ref['clot_count']);ids=np.arange(len(X)) if order is None else order
 model=gs.Gaussian(dim=3,var=1.,len_scale=LC,rescale=1.,nugget=0.)
 srf=gs.SRF(model,seed=SEED,generator='RandMeth',mode_no=4096)
 z=np.empty(len(X));z[ids]=srf(X[ids].T,mesh_type='unstructured')
 shift=brentq(lambda c:np.average(np.tanh(z[:n]-c),weights=ref['V'][:n]),z.min()-10,z.max()+10,xtol=1e-14)
 q=np.tanh(z-shift);b=ref['bonds'];bulk=ref['bond_types']==1;qb=q[b].mean(1)
 w=ref['taper']*ref['r0']*V**2*(1/ref['m'][b[:,0]]+1/ref['m'][b[:,1]])
 Z=np.average(np.exp(ALPHA*qb[bulk]),weights=w[bulk]);sc=SC[ref['bond_types']].copy()
 sc[bulk]=S0[1]+(SC[1]-S0[1])*np.exp(ALPHA*qb[bulk])/Z
 return dict(atom_ids=np.arange(1,len(X)+1,dtype=np.int64),reference_X=X,q=q,z=z,center_shift=np.array(shift),bond_ids=np.arange(1,len(b)+1,dtype=np.int64),bond_atom_ids=b+1,bond_types=ref['bond_types'],bond_q=qb,bond_sc=sc,bond_s0=S0[ref['bond_types']],bond_weight=w,Z=np.array(Z),reference_volume=ref['V'],node_bulk_sc_proxy=S0[1]+(SC[1]-S0[1])*np.exp(ALPHA*q)/Z)
def main():
 ref=load_npz(O/'provenance/reference.npz');a=generate(ref);again=generate(ref,np.arange(len(ref['X']))[::-1]);assert all(np.array_equal(a[k],again[k]) for k in a),'ID-order reproducibility failed'
 assert not (O/'HETEROGENEITY_FIELD.npz').exists();np.savez_compressed(O/'HETEROGENEITY_FIELD.npz',**a);dump(O/'HETEROGENEITY_SPEC.json',SPEC)
 n=int(ref['clot_count']);b=ref['bonds'];bulk=ref['bond_types']==1;w=a['bond_weight'];q=a['q'][:n];X=ref['X'][:n]
 axes=[np.unique(X[:,c]) for c in range(3)];idx=tuple(np.searchsorted(axes[c],X[:,c]) for c in range(3));grid=np.empty(tuple(map(len,axes)));grid[idx]=q
 interior=np.ones(grid.shape,bool)
 for c in range(3):
  sl=[slice(None)]*3;sl[c]=0;interior[tuple(sl)]=False;sl[c]=-1;interior[tuple(sl)]=False
 extrema={}
 for kind,flt in [('maxima',maximum_filter),('minima',minimum_filter)]:
  mask=(grid==flt(grid,size=3,mode='nearest'));extrema[kind]=dict(all_nodes=int(mask.sum()),interior_nodes=int((mask&interior).sum()))
 corr=[]
 for c in range(3):
  for lag in range(len(axes[c])):
   sl1=[slice(None)]*3;sl2=[slice(None)]*3
   if lag:sl1[c]=slice(0,-lag);sl2[c]=slice(lag,None)
   x=grid[tuple(sl1)].ravel();y=grid[tuple(sl2)].ravel();corr.append(dict(axis='XYZ'[c],distance_m=lag*DX,empirical_q_covariance_normalized=float(np.mean(x*y)/np.mean(q*q)),pairs=len(x),nominal_z_covariance=float(np.exp(-(lag*DX/LC)**2))))
 pred=np.c_[np.ones(n),X]@np.linalg.lstsq(np.c_[np.ones(n),X],q,rcond=None)[0];r2=1-np.sum((q-pred)**2)/np.sum((q-q.mean())**2)
 meta=dict(**SPEC,gstools_version=gs.__version__,field_sha256=sha(O/'HETEROGENEITY_FIELD.npz'),Z=float(a['Z']),q_min=float(q.min()),q_max=float(q.max()),q_volume_mean=float(np.average(q,weights=ref['V'][:n])),q_std=float(q.std()),linear_gradient_R2=float(r2),extrema_26_neighbor=extrema,
    axis_span_m=np.ptp(X,axis=0),lc_over_dx=LC/DX,lc_over_horizon=LC/HORIZON,horizon_over_dx=HORIZON/DX,
    axis_slice_mean_ranges={c:[float(grid.mean(axis=tuple(k for k in range(3) if k!=i)).min()),float(grid.mean(axis=tuple(k for k in range(3) if k!=i)).max())] for i,c in enumerate('XYZ')},
    sc_min=float(a['bond_sc'][bulk].min()),sc_max=float(a['bond_sc'][bulk].max()),sc_weighted_mean=float(np.average(a['bond_sc'][bulk],weights=w[bulk])),sc_uncut_weighted_mean=float(np.average(a['bond_sc'][bulk&~ref['initial_cut']],weights=w[bulk&~ref['initial_cut']])),
    positive_q_components_26=int(label(grid>0,np.ones((3,3,3)))[1]),negative_q_components_26=int(label(grid<0,np.ones((3,3,3)))[1]),
    correlation_warning='Nominal Gaussian z covariance is not empirical finite-domain transformed q covariance; extrema and excursion domains are not statistically independent tissue regions.',reproducible_reverse_ID_order=True)
 dump(O/'HETEROGENEITY_METADATA.json',meta)
 for path,rows in [('FIELD_CORRELATION.csv',corr),('NODE_LOCAL_PROPERTIES.csv',[dict(atom_id=i+1,role='material' if i<n else 'anchor',X_m=a['reference_X'][i,0],Y_m=a['reference_X'][i,1],Z_m=a['reference_X'][i,2],q=a['q'][i],local_bulk_sc_proxy=a['node_bulk_sc_proxy'][i],reference_volume_m3=ref['V'][i]) for i in range(len(a['q']))]),('BOND_LOCAL_PROPERTIES.csv',[dict(original_bond_id=i+1,atom_i=int(b[i,0]+1),atom_j=int(b[i,1]+1),type=int(ref['bond_types'][i]),initial_precut=bool(ref['initial_cut'][i]),q_ij=a['bond_q'][i],s0=a['bond_s0'][i],sc=a['bond_sc'][i],Delta_s=a['bond_sc'][i]-a['bond_s0'][i],reference_weight=w[i]) for i in range(len(b))])]:
  with (O/path).open('w') as f:wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
 print('FIELD_FROZEN',meta,flush=True)
if __name__=='__main__':main()
