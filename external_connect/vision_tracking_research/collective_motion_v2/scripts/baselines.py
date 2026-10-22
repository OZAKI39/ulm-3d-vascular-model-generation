"""Re-run frozen F10U libraries on the identical F10V raw endpoint images."""
import sys,argparse,time,json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,cv2
from legacy import piv,cpd,CandidateDetector
from vcommon import OUT,DATA,P,NAMES,dump,identity,stats,sha
A=None;DT=None;DET=None;METHOD=None

def init(split,method):
 global A,DT,DET,METHOD
 cv2.setNumThreads(1);A=np.load(DATA/split/'raw.npy',mmap_mode='r');DT=np.load(DATA/split/'dt.npy',mmap_mode='r');DET=CandidateDetector();METHOD=method

def worker(i):
 t=time.perf_counter();oa,_=DET(A[i,0]);ob,_=DET(A[i,2]);detect=(time.perf_counter()-t)*1000
 if METHOD=='CPD':r=cpd([o['xy'] for o in oa],[o['xy'] for o in ob],sum(DT[i]),A[i,0].shape)
 else:r=piv(A[i,0],A[i,2],float(sum(DT[i])),mode='candidate',objects_a=oa,objects_b=ob,multipass=METHOD=='PIV_MULTIPASS')
 e=r['estimate'];mean=np.array(e.get('raw_estimate',[np.nan]*6));accepted=all(e['valid'].values());reason=next((v for v in e['reasons'].values() if v is not None),None)
 label=0 if accepted else (1 if reason=='NON_AFFINE_OR_INCONSISTENT_VECTORS' else 2)
 residual=e.get('residual');score=float(residual) if residual is not None else 1e6
 row=dict(index=i,method=METHOD,mean=mean,accepted=accepted,score=score,predicted_class=label,estimate=e,detect_endpoints_ms=detect,pipeline_ms=(time.perf_counter()-t)*1000,vectors=r.get('vectors',[]),algorithm_ms=r['algorithm_ms'])
 return row

def run(split,method,workers):
 if split=='final_test' and not (OUT/'results/model_freeze.json').exists():raise RuntimeError('FINAL_TEST_LOCKED')
 n=len(np.load(DATA/split/'dt.npy',mmap_mode='r'));indices=list(range(n)) if method=='PIV_SINGLE' else list(range(min(112,n)))
 start=time.perf_counter()
 with ProcessPoolExecutor(max_workers=workers,initializer=init,initargs=(split,method)) as pool:rows=list(pool.map(worker,indices,chunksize=8))
 file=OUT/'results'/f'{method}_{split}.json';dump(file,dict(rows=rows,identity=identity(),dataset_manifest_sha256=sha(DATA/split/'manifest.json'),wall_s=time.perf_counter()-start,parallel_workers=workers,scope='algorithm accuracy; per-call timings under CPU concurrency, use separate runtime benchmark for fair latency',pipeline_ms=stats([x['pipeline_ms'] for x in rows])))
 print(method,split,len(rows),'wall',time.perf_counter()-start,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--split',choices=['validation','final_test'],required=True);p.add_argument('--methods',nargs='+',default=['PIV_SINGLE','PIV_MULTIPASS','CPD']);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
 for m in a.methods:run(a.split,m,a.workers)
