import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *
from affine import parameters
import h5py

p=SIM/'results/straight_tube_exp1/f10s/a_navigation_v1/runs/X_PAIR_3MS_DT1/trajectory.h5'
gate=json.loads((SIM/'results/straight_tube_exp1/f10s/analysis/FINAL_GATE.json').read_text())['groups']['X_PAIR_3MS_DT1']
assert gate['simulation_group_position_eligible'] and gate['simulation_group_shape_eligible']
with h5py.File(p) as f:t=f['time_s'][:];pos=f['positions_m'][:];ids=f['particle_id'][:]
xy=pos[:,:,:2]*1e6; x=xy[0];y=xy[-1];center=x.mean(0);X=np.c_[np.ones(len(x)),x-center];B=np.linalg.lstsq(X,y-x,rcond=None)[0]
F=np.eye(2)+B[1:].T;b=B[0]-(F-np.eye(2))@center
par=parameters(F,b,t[-1]-t[0],center);res=y-x-X@B
regions={}
for name,mask in [('upper',x[:,1]>=center[1]),('lower',x[:,1]<center[1])]:
    if mask.sum()>=4:
        BB=np.linalg.lstsq(X[mask],(y-x)[mask],rcond=None)[0]
        regions[name]={'n':int(mask.sum()),'affine_displacement_coeff_um':BB.tolist(),'residual_rms_um':float(np.sqrt(np.mean(((y-x)[mask]-X[mask]@BB)**2)))}
rows={'identity':identity(),'file':str(p),'sha256':sha(p),'source_git':'31a44d3f342e5d0bca55f7f9b83f4d32b8189696',
      'time_samples':len(t),'particle_count':len(ids),'window_s':float(t[-1]-t[0]),'velocity_source':'finite differences of saved positions / saved time; velocities not stored as dataset',
      'com_displacement_um':(pos[-1].mean(0)-pos[0].mean(0))*1e6,'xy_covariance_start_um2':np.cov(x.T,bias=True),'xy_covariance_end_um2':np.cov(y.T,bias=True),
      'six_parameters':dict(zip(NAMES,par)),'translation_unit':'um/s','gradient_unit':'1/s','global_residual_rms_um':np.sqrt(np.mean(res**2)),
      'local_regions':regions,'gate':gate,'interpretation':'Fixed IDs used ONLY as saved simulator ground truth; no vision identities. Projection loses out-of-plane motion. Small 3 ms shape signal is not evidence of arbitrary deformation control.',
      'synthetic_image_rendering':'NOT_RUN; no microscope PSF or pixel calibration for this device',
      'PHYSICS_UNCALIBRATED':True,'training_ready':False,'new_simulation_run':False}
dump(ROOT/'results/simulation_bridge.json',rows)
np.savez_compressed(ROOT/'results/simulation_projection.npz',time_s=t,xy_um=xy,velocity_um_s=np.diff(xy,axis=0)/np.diff(t)[:,None,None])
print('simulation projection:',rows['com_displacement_um'],par)
