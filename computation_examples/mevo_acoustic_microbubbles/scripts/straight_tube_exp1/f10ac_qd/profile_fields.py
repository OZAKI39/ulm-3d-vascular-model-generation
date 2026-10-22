"""Short relevant field-call profile, explicitly synchronized transfer/compute phases."""
import time,h5py,numpy as np,warp as wp
from scripts.straight_tube_exp1.f10ab.run import AuditEngine
from scripts.straight_tube_exp1.f10ad_multifrequency.interface import context,make_schedule
from .common import *
def main():
 c=context()
 with h5py.File(BASE/'INITIAL_STATES.h5') as h:x=h['S2_r0'][:]
 e=AuditEngine(x,c['physics'],c['response'],c['basis'],c['cavity'],policy=policy());e.control(make_schedule([0,0,700,700],[0,0,0,np.pi]).segments[0].event());e.fields(x);s=e.sim;rows=[]
 for _ in range(12):
  ts=[time.perf_counter()];s.x.assign(x);wp.synchronize_device(s.device);ts.append(time.perf_counter())
  s.fields();wp.synchronize_device(s.device);ts.append(time.perf_counter())
  p=s.pressure.numpy();f1=s.force.numpy().copy();ts.append(time.perf_counter())
  s.pairs();wp.synchronize_device(s.device);ts.append(time.perf_counter())
  f=s.force.numpy().copy();bits=int(s.error.numpy()[0]);ts.append(time.perf_counter());assert bits==0
  rows.append(np.diff(ts))
 write(BASE/'evidence/FIELD_TRANSFER_PROFILE.json',dict(columns=['position_host_to_device','primary_kernel_and_launch','pressure_and_primary_device_to_host','secondary_kernel_and_launch','total_and_error_device_to_host'],samples_s=rows,mean_s=np.mean(rows,axis=0),N=len(x),note='12 useful fixed-state profiles after warmup; synchronization added only to this diagnostic, no trajectory code change; not total end-to-end step cost',identity=identity()))
if __name__=='__main__':main()
