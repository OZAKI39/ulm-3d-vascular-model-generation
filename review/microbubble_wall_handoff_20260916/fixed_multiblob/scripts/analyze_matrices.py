"""Independent analysis: never imports the hydrodynamic solver adapter."""
from pathlib import Path
import json, numpy as np,h5py
R=Path(__file__).resolve().parents[1]
C=json.loads((R/'contracts/FIXED_MULTIBLOB_WALL_FEASIBILITY_CONTRACT.json').read_text())
with h5py.File(R/'input/RMBW_WALL_RESISTANCE_TABLE_V0.h5') as f:EPS=f['epsilon'][:];REF=f['R_total_scaled'][:]
Q=np.random.default_rng(C['seeds']['actions']).normal(size=(256,6));Q/=np.linalg.norm(Q,axis=1)[:,None]
D=np.random.default_rng(C['seeds']['dissipation']).normal(size=(10000,6));D/=np.linalg.norm(D,axis=1)[:,None]
def reference(e):return np.array([np.interp(np.log(e),np.log(EPS),REF[:,i,j]) for i in range(6) for j in range(6)]).reshape(6,6)
def modes(A):return dict(normal=A[2,2],tangent=(A[0,0]+A[1,1])/2,RR_parallel=(A[3,3]+A[4,4])/2,RR_normal=A[5,5],TR=(A[1,3]-A[0,4])/2)
def props(A):
 ev=np.linalg.eigvalsh((A+A.T)/2);return dict(reciprocity=float(np.linalg.norm(A-A.T)/np.linalg.norm(A)),min_eigenvalue=float(ev.min()),condition=float(np.linalg.cond(A)),minimum_dissipation=float(np.min(np.einsum('ni,ij,nj->n',D,A,D))))
def compare(A,e):
 B=reference(e);ma=modes(A);mb=modes(B);out={k:float(v) for k,v in ma.items()}
 for k in ma:out[k+'_error']=float(abs(ma[k]-mb[k])/max(abs(mb[k]),1e-12))
 out['matrix_action_error']=float(np.max(np.linalg.norm(Q@(A-B).T,axis=1)/np.linalg.norm(Q@B.T,axis=1)))
 out['excess_matrix_action_error']=float(np.max(np.linalg.norm(Q@(A-B).T,axis=1)/np.linalg.norm(Q@(B-np.eye(6)).T,axis=1)))
 out['TT_normal_error']=float(abs(A[2,2]-B[2,2])/B[2,2]);out['TT_parallel_error']=float(max(abs(A[0,0]-B[0,0]),abs(A[1,1]-B[1,1]))/B[0,0]);out['RR_error']=max(out['RR_parallel_error'],out['RR_normal_error']);return out

def extract(paths):
 rows=[];geometry=[];leaks=[];performance=[];failures=[]
 for fn in paths:
  with h5py.File(fn) as f:
   for gn,g in f.items():
    meta={k:(v.item() if isinstance(v,np.generic) else v) for k,v in g.attrs.items()};meta['geometry']=gn;meta['raw_file']=fn.name;meta['file_terminal_state']=str(f.attrs.get('terminal_state','NO_HDF_TERMINAL_RECORD_CHECK_EXECUTION_RECEIPT'));geometry.append(meta)
    for qn,q in g.items():
     if not isinstance(q,h5py.Group):continue
     if q.attrs.get('status')!='COMPLETED':failures.append(dict(**meta,query=qn,reason=q.attrs.get('reason','unknown')));continue
     if 'epsilon' in q.attrs:
      row={k:meta[k] for k in ['geometry','raw_file','category','layout','layers','extent','spacing_ratio','N_wall','beta']};row.update({k:(v.item() if isinstance(v,np.generic) else v) for k,v in q.attrs.items()});row.update(query=qn,x=float(q['targets'][0,0]),y=float(q['targets'][0,1]),z=float(q['targets'][0,2]));A=q['R_scaled'][:];row.update(compare(A,row['epsilon']),**props(A));rows.append(row)
    if 'LEAKAGE' in g and 'M_scaled' in g['LEAKAGE'] and 'M_scaled' in g['LEAKAGE_FREE']:
     M=g['LEAKAGE/M_scaled'][:];M0=g['LEAKAGE_FREE/M_scaled'][:];leaks.append(dict(**meta,transmission=float(M[8,2]/M0[8,2]),upper_force=1.,lower_velocity_scaled=float(M[8,2]),free_lower_velocity_scaled=float(M0[8,2])))
    if meta['category']=='PERFORMANCE':
     vals=[q for q in g.values() if isinstance(q,h5py.Group) and 'R_scaled' in q];performance.append(dict(**meta,query_count=len(vals),median_query_seconds=float(np.median([q.attrs['total_seconds'] for q in vals])) if vals else None,median_assembly_seconds=float(np.median([q.attrs['assembly_seconds'] for q in vals])) if vals else None,median_solve_seconds=float(np.median([q.attrs['solve_seconds'] for q in vals])) if vals else None,max_rss_bytes=max([q.attrs.get('peak_rss_bytes',0) for q in g.values() if isinstance(q,h5py.Group)],default=None),first_including_precompute_seconds=float(meta.get('wall_total_seconds',0)+vals[0].attrs['total_seconds']) if vals else None))
 return rows,geometry,leaks,performance,failures

def groups(rows,keys):
 out={}
 for r in rows:out.setdefault(tuple(r[k] for k in keys),[]).append(r)
 return out

def phase_statistics(rows):
 out=[]
 keys=['beta','layout','layers','extent','spacing_ratio','epsilon']
 for key,rr in groups([r for r in rows if r['category']=='PLANAR'],keys).items():
  for m in ['normal','tangent','RR_parallel','RR_normal','TR']:
   v=np.array([r[m] for r in rr]);mean=float(v.mean());out.append(dict(zip(keys,key),mode=m,mean=mean,std=float(v.std()),minimum=float(v.min()),maximum=float(v.max()),CV=float(v.std()/max(abs(mean),1e-12)),spread=float(np.ptp(v)/max(abs(mean),1e-12)),phase_count=len(v)))
 return out
if __name__=='__main__':
 import sys
 rows,geo,leaks,perf,fail=extract([Path(p) for p in sys.argv[1:]]);phase=phase_statistics(rows)
 for b in sorted(set(r['beta'] for r in rows)):
  core=[r for r in rows if r['beta']==b and .01<=r['epsilon']<=.2 and r['category']=='PLANAR'];ps=[r for r in phase if r['beta']==b and .01<=r['epsilon']<=.2 and r['mode'] in ['normal','tangent']]
  print(json.dumps(dict(beta=b,queries=len(core),worst_action=max(r['matrix_action_error'] for r in core),mean_action=float(np.mean([r['matrix_action_error'] for r in core])),worst_TT_phase=max(r['spread'] for r in ps)),indent=2))
 print('FAILURES',len(fail))
