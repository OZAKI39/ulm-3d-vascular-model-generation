"""Match exact HF nodes to recorded LF nodes/dense samples, without HF interpolation."""
import numpy as np,h5py
from .common import *
def aligned(h,z):
 ht=h['time_s'][:];candidates=[]
 for tk,xk in [('sample_time_s','sample_positions_m'),('time_s','positions_m')]:
  lt=z[tk][:];ix=np.searchsorted(ht,lt);ix=np.clip(ix,0,len(ht)-1);left=np.maximum(ix-1,0);ix=np.where(abs(ht[left]-lt)<abs(ht[ix]-lt),left,ix);good=np.flatnonzero(abs(ht[ix]-lt)<1e-12);_,unique=np.unique(ix[good],return_index=True);good=good[unique]
  candidates.append((len(good),tk,xk,good,ix))
  if len(good)>=101:break
 count,tk,xk,good,ix=max(candidates,key=lambda c:c[0]);assert count>=min(101,len(ht)),'Insufficient exact shared times; no numeric admission'
 good=good[np.unique(np.rint(np.linspace(0,len(good)-1,min(401,len(good)))).astype(int))];ids=ix[good];return ht[ids],h['positions_m'][ids],z[xk][good],dict(LF_source=tk,shared_samples=len(good),HF_interpolation=False,scope='Exact shared timestamps; original LF dense output or accepted nodes')
def compare(t):
 a=read(BASE/'hf'/t['id']/'SUMMARY.json');b=read(BASE/'lf'/t['id']/'SUMMARY.json');r=dict(id=t['id'],HF_status=a['status'],LF_status=b['status'],PASS=False)
 if a['status']!='PASS' or b['status']!='PASS':return r
 with h5py.File(BASE/a['trajectory_path']) as h,h5py.File(BASE/b['trajectory_path']) as z:
  tt,xx,yy,meta=aligned(h,z);err=np.linalg.norm(xx-yy,axis=2)
 r.update(max_fixed_ID_m=float(err.max()),PASS=bool(err.max()<=50e-9),paired_duration_s=float(tt[-1]),alignment=meta);return clean(r)
