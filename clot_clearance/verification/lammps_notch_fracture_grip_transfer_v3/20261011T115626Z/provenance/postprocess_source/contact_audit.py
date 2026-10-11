"""Read-only source-equivalent conservative pair contact on saved accepted states."""
from pathlib import Path
import sys,math,numpy as np
from scipy.spatial import cKDTree
from design import save,js
from analyze import loadcsv

def main(base):
 p=base/'M5C';g=dict(np.load(p/'geometry.npz'));b=g['bonds'];X=g['X'];N=len(X);ev=np.load(p/'break_events.npz')['events'];bt=np.full(len(b),2**62,np.int64);bt[ev[:,1]]=ev[:,0];a=loadcsv(p/'history.csv');V=(50e-6)**3;delta=150.5e-6;spring=15*(18*(100000/3)/(math.pi*delta**4))*V*V/delta;records=[]
 for path in sorted(p.glob('*.npz')):
  z=dict(np.load(path))
  if 'x_m' not in z:continue
  x=z['x_m'];step=int(z['step']);pairs=cKDTree(x).query_pairs(50e-6,output_type='ndarray');key=pairs[:,0]*N+pairs[:,1];oldkey=b[bt>=step,0]*N+b[bt>=step,1];pairs=pairs[~np.isin(key,oldkey)];d=x[pairs[:,0]]-x[pairs[:,1]];r=np.linalg.norm(d,axis=1);gap=50e-6-r;force=spring*gap[:,None]*d/r[:,None] if len(pairs) else np.empty((0,3));f=np.zeros_like(x);np.add.at(f,pairs[:,0],force);np.add.at(f,pairs[:,1],-force);U=.5*spring*np.sum(gap*gap);native=float(a['U_contact_J'][step]);err=abs(U-native);rec=dict(file=path.name,step=step,time_s=step*.5e-6,active_pairs=len(pairs),max_node_force_N=float(np.linalg.norm(f,axis=1).max()),max_pair_force_N=float(np.linalg.norm(force,axis=1).max()) if len(pairs) else 0,source_equivalent_U_J=float(U),native_U_J=native,absolute_U_error_J=float(err),relative_U_error=float(err/max(native,1e-30)),top_contact_Rz_N=float(-f[g['prescribed_mask'][:,2]&(X[:,2]>0),2].sum()),bottom_contact_Rz_N=float(-f[g['prescribed_mask'][:,2]&(X[:,2]<0),2].sum()))
  records.append(rec)
  if path.stem=='end':np.savez_compressed(base/'analysis/M5C/end_contact_fields.npz',pairs=pairs,force_pairs_N=force,force_nodes_N=f,x_m=x,step=step)
 passed=all(q['absolute_U_error_J']<1e-22+1e-8*q['native_U_J'] for q in records)
 result=dict(status='SOURCE_LAW_AND_NATIVE_ENERGY_AGREE_AT_SAVED_STATES' if passed else 'NOT_VALIDATED_DO_NOT_USE_FORCE_AS_NATIVE_MEASUREMENT',source='Locked pair_bpm_peri.cpp: F=15 k V^2 (dx-r)/delta; r<dx; exclude currently bonded pairs before this step bond deletion; no friction.',scope='Reconstruction on every retained native snapshot only, not an every-step maximum. Native pair compute precedes bond rupture; exclusion uses previous accepted alive set. Energy comparison tests the actual candidate/exclusion set.',records=records,max_force_over_saved_states_N=max(q['max_node_force_N'] for q in records),max_native_energy_error_J=max(q['absolute_U_error_J'] for q in records));save(base/'CONTACT_FORCE_AUDIT.json',result);print(result)
 assert passed,'Contact force reconstruction did not match native energy; retain explicit unavailable status'

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
