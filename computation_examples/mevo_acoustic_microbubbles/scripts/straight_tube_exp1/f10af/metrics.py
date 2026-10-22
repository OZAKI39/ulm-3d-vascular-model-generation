"""Fixed-ID 3D truth and the unmodified F10Y 2D distribution contract."""
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import cdist
from .vendor.f10y_aggregation import distribution
from scripts.straight_tube_exp1.f10ab.metrics import measure
from .common import R

def spatial(x):
    x=np.asarray(x);q=x[:,:2];extent=np.ptp(q,axis=0);theta=np.radians(np.arange(0,180,5));freq=np.linspace(2000,6500,55)
    k=(2*np.pi*freq[:,None,None]*np.c_[np.cos(theta),np.sin(theta)][None,:,:]).reshape(-1,2)
    # Fixed 1.0 x 0.8 mm observation window. Remove known rectangular support response.
    window=np.array([.001,.0008]);phase=q@k.T;null=np.prod(np.sinc(k*window/(2*np.pi)),axis=1);amp=np.exp(1j*phase).mean(0);S=len(x)*abs(amp-null)**2
    order=np.argsort(S)[::-1];j=order[0];vec=k[j];spacing=2*np.pi/np.linalg.norm(vec);direction=vec/np.linalg.norm(vec);u=q@direction;periods=float(np.ptp(u)/spacing);bands=np.rint(u/spacing-np.angle(amp[j])/(2*np.pi)).astype(int);counts=np.unique(bands,return_counts=True)[1]
    occupied=int(np.count_nonzero(counts>=max(2,int(.03*len(x)))))
    peak=float(S[j]);stripe=peak>=12 and periods>=2 and occupied>=3
    other=[i for i in order if abs(np.dot(k[i]/np.linalg.norm(k[i]),direction))<.5]
    second=other[0];periods2=float(np.ptp(q@(k[second]/np.linalg.norm(k[second])))*np.linalg.norm(k[second])/(2*np.pi))
    lattice=bool(stripe and S[second]>=12 and periods2>=2)
    centered=q-q.mean(0);rad=np.linalg.norm(centered,axis=1);angle=np.arctan2(centered[:,1],centered[:,0]);coverage=np.count_nonzero(np.histogram(angle,bins=np.linspace(-np.pi,np.pi,13))[0])/12
    ring_cv=float(rad.std()/max(rad.mean(),1e-30));center_hole=float(rad.min()/max(rad.mean(),1e-30));nn=cKDTree(q).query(q,k=min(7,len(x)));delta=q[nn[1][:,1]]-q;angles=np.arctan2(delta[:,1],delta[:,0]);psi4=abs(np.exp(4j*angles).mean());psi6=abs(np.exp(6j*angles).mean())
    c=x.mean(0);cov=(x-c).T@(x-c)/len(x);comp,labels=connected_components(cdist(x,x)<=2.1*R,directed=False);localcomp,locallabel=connected_components(cdist(x,x)<=20e-6,directed=False)
    xy=distribution(q*1e6)
    return dict(N=len(x),center_m=c,Rg_3d_m=float(np.sqrt(np.trace(cov))),covariance_3d_m2=cov,xy_f10y=xy,xy_coordinate_unit='micrometre numerical coordinates, not calibrated camera pixels',extent_m=np.ptp(x,axis=0),minimum_gap_m=float(cdist(x,x)[np.triu_indices(len(x),1)].min()-2*R),proximity_components=int(comp),largest_component=int(np.bincount(labels).max()),local_components_20um=int(localcomp),largest_local_component=int(np.bincount(locallabel).max()),structure_peak_S=peak,stripe_order=float(abs(amp[j])),stripe_spacing_m=float(spacing),stripe_normal_deg=float(np.degrees(np.arctan2(vec[1],vec[0]))%180),periods=periods,occupied_bands=occupied,stripe_descriptor_pass=bool(stripe),second_peak_S=float(S[second]),lattice_descriptor_pass=lattice,psi4=float(psi4),psi6=float(psi6),ring_cv=ring_cv,center_hole=center_hole,angular_coverage=coverage,ring_descriptor_pass=bool(ring_cv<.15 and center_hole>.5 and coverage>=.75),spectrum_S=S.reshape(len(freq),len(theta)),spectrum_freq_m1=freq,spectrum_theta_rad=theta)
def rotation(x,v):
    r=x-x.mean(0);w=v-v.mean(0);den=np.sum(r*r);denxy=np.sum(r[:,:2]**2);cross=np.cross(r,w)[:,2].sum();omega=cross/max(den,1e-99);omxy=cross/max(denxy,1e-99)
    fit=np.linalg.lstsq(np.c_[np.ones(len(x)),r[:,:2]],w[:,:2],rcond=None)[0];A=fit[1:].T;rigid=np.c_[-omxy*r[:,1],omxy*r[:,0],np.zeros(len(r))];energy=np.sum(w*w);fraction=1-np.sum((w-rigid)**2)/max(energy,1e-99) if energy>1e-40 else 0.
    return dict(Omega_z_s=float(omega),Omega_xy_s=float(omxy),F10Y_omega_affine_s=float(A[1,0]-A[0,1]),affine_A_s=A,shear_s=float(np.linalg.norm((A+A.T)/2-np.trace(A)*np.eye(2)/2)),rigid_fraction=float(fraction),mean_velocity_m_s=v.mean(0),relative_velocity_rms_m_s=float(np.sqrt(energy/len(x))),rotation_center_m=x.mean(0),rotation_center_scope='Group COM reference; no inferred fixed physical vortex center')
def endpoint(x0,x,v=None):
    out=spatial(x);m=measure(x0,x);out.update(COM_displacement_m=float(np.linalg.norm(x.mean(0)-x0.mean(0))),max_fixed_id_displacement_m=float(np.linalg.norm(x-x0,axis=1).max()),relative_fixed_ID_rms_m=float(m['translation_removed_rms_3d_m']),kabsch_rms_m=float(m['kabsch_rmsd_3d_m']),pair_distance_rms_change_m=float(m['pair_distance_rms_change_3d_m']),rigid_rotation_3d=m['rigid_rotation_3d'])
    if v is not None:out['rotation']=rotation(x,v)
    return out
