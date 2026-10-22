"""Official CuPy tutorial plus known translation, blank, and backend identity checks."""
import sys,time,traceback,json,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from backends import Backend
import cv2
def main():
 cv2.setNumThreads(2);rows=[];fail=[];mon=ResourceMonitor()
 adapters={n:Backend(n) for n in ['B0','B1','B2','B3','B4']}
 rng=np.random.default_rng(201017);a=(cv2.GaussianBlur(rng.uniform(0,255,(240,240)).astype(np.float32),(0,0),1.2)).astype(np.uint8)
 for case,dx,dy in [('static',0,0),('x_positive',2,0),('y_negative',0,-3),('subpixel',.4,.7),('large',10,-8),('blank',0,0)]:
  aa=np.full_like(a,180) if case=='blank' else a
  bb=cv2.warpAffine(aa,np.array([[1,0,dx],[0,1,dy]],np.float32),(240,240),borderMode=cv2.BORDER_REFLECT)
  outputs={}
  for name,engine in adapters.items():
   try:
    r=engine.run(aa,bb,.0666,mode='raw');outputs[name]=r;good=r['valid'];e=np.hypot(r['u']-dx,r['v']-dy)
    rows.append(dict(case=case,name=name,truth_dxdy=[dx,dy],accepted_vectors=int(good.sum()),n_vectors=good.size,mae_epe_all=float(np.nanmean(e)),mae_epe_valid=float(np.mean(e[good])) if good.any() else None,estimate=r['estimate'],dtype=r['dtype']))
   except Exception as e:fail.append(dict(case=case,name=name,error=repr(e),traceback=traceback.format_exc()))
  if 'B1' in outputs and 'B2' in outputs:
   r,s=outputs['B1'],outputs['B2'];joint=np.isfinite(r['u']+r['v']+s['u']+s['v'])
   rows.append(dict(case=case,name='RUST_VS_SCIPY',max_component_difference=float(np.max(np.abs(np.stack([r['u']-s['u'],r['v']-s['v']]))[:,joint])) if joint.any() else None,valid_disagreements=int(np.count_nonzero(r['valid']!=s['valid'])),snr_max_difference=float(np.nanmax(np.abs(r['snr']-s['snr']))) if joint.any() else None))
 # Reproduce the official basic notebook's first case and default multipass.
 engine=adapters['B4'];folder=THIRD/'openpiv-python-gpu/openpiv_gpu/tutorials/test0';paths=sorted(folder.glob('*.tiff'))
 try:
  aa,bb=[cv2.imread(str(p),cv2.IMREAD_ANYDEPTH) for p in paths];t=time.perf_counter()
  solver=engine.gp.piv_gpu(aa.shape,min_search_size=16,search_size_iters=(1,2,2));u,v=solver(aa,bb);engine.sync()
  official=dict(status='COMPLETED',wall_s=time.perf_counter()-t,input_files={str(p):sha(p) for p in paths},shape=aa.shape,field_shape=u.shape,n_finite=int(np.isfinite(u+v).sum()),replaced_or_rejected=int(solver.val_locations.sum()),note='Official tutorial numerical execution, including native replacement. Demonstration only; replacements are not admitted measurements.')
  np.savez_compressed(OUT/'results/cupy_official_example.npz',u=u,v=v,x=solver.coords[0],y=solver.coords[1],replaced=solver.val_locations)
 except Exception as e:official=dict(status='FAILED',error=repr(e),traceback=traceback.format_exc())
 sources=[]
 for p in THIRD.iterdir():
  if (p/'.git').exists():sources.append(dict(name=p.name,head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=p,text=True).strip(),diff=subprocess.check_output(['git','diff'],cwd=p,text=True),licenses={f.name:sha(f) for f in p.iterdir() if 'LICEN' in f.name.upper()}))
 dump(OUT/'results/official_and_consistency.json',dict(identity=identity(),rows=rows,failures=fail,cupy_official=official,sources=sources,resources=mon.close()))
 print(json.dumps(clean(dict(rows=rows,failures=fail,cupy_official=official))),flush=True)
if __name__=='__main__':main()
