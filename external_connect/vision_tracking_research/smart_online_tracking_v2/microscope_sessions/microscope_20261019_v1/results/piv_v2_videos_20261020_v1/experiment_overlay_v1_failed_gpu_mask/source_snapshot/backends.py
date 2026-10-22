"""Thin adapters only: all correlation, peak fitting and validation are upstream."""
import importlib.util,inspect,time,sys,contextlib,warnings
import numpy as np
import cv2
from xcommon import THIRD,fit_affine
from openpiv import pyprocess,validation

def preprocess(a,mode='local',background=None):
 if mode=='raw':return np.asarray(a,np.float32)
 if mode=='local':return cv2.subtract(cv2.GaussianBlur(a,(0,0),4),a).astype(np.float32)
 if mode=='clahe':return cv2.createCLAHE(clipLimit=2.,tileGridSize=(8,8)).apply(a).astype(np.float32)
 if mode in ['fixed','causal']:
  if background is None:raise ValueError('BACKGROUND_REQUIRED')
  return np.maximum(np.asarray(background,np.float32)-a,0)
 raise ValueError(mode)

def window_std(a,x,y,ws):
 # OpenCV integral images replace Python per-window reductions, not PIV kernels.
 ii,ss=cv2.integral2(np.asarray(a,np.float32),sdepth=cv2.CV_64F,sqdepth=cv2.CV_64F)
 x0=np.clip(np.rint(x-ws/2).astype(int),0,a.shape[1]-ws);y0=np.clip(np.rint(y-ws/2).astype(int),0,a.shape[0]-ws)
 x1=x0+ws;y1=y0+ws
 mean=(ii[y1,x1]-ii[y0,x1]-ii[y1,x0]+ii[y0,x0])/(ws*ws)
 var=(ss[y1,x1]-ss[y0,x1]-ss[y1,x0]+ss[y0,x0])/(ws*ws)-mean**2
 return np.sqrt(np.maximum(var,0))

class Backend:
 def __init__(self,name,ws=48,ov=24,profile=False):
  self.name=name;self.ws=ws;self.ov=ov;self.cache={};self.profile=profile;self.parts={}
  self.pp=pyprocess
  if name=='B0':
   path='/workspace/f10w_vision_20261017/.venv/lib/python3.12/site-packages/openpiv/pyprocess.py'
   spec=importlib.util.spec_from_file_location('f10w_frozen_pyprocess',path);self.pp=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.pp)
  if name=='B2':
   import openpiv_rust
   self.rust_path=openpiv_rust.__file__
  if name=='B3':
   sys.path.insert(0,str(THIRD/'TorchPIV/src'))
   import torch
   # Upstream package __init__ eagerly imports its GUI and forces Qt5Agg.
   # Load the unchanged public compute module through a namespace package:
   # no Qt window and no replacement of any scientific function.
   import types
   package=types.ModuleType('torchPIV');package.__path__=[str(THIRD/'TorchPIV/src/torchPIV')]
   sys.modules.setdefault('torchPIV',package)
   from torchPIV import PIVbackend
   torch.set_num_threads(2);self.torch=torch;self.tp=PIVbackend
  if name=='B4':
   sys.path.insert(0,str(THIRD/'openpiv-python-gpu'))
   import cupy as cp
   from openpiv_gpu import gpu_process
   self.cp=cp;self.gp=gpu_process
 def sync(self):
  if self.name=='B3':self.torch.cuda.synchronize()
  if self.name=='B4':self.cp.cuda.get_current_stream().synchronize()
 @contextlib.contextmanager
 def instrument(self):
  self.parts={};original=[]
  if not self.profile:yield;return
  targets=[]
  if self.name in ['B0','B1','B2']:targets=[(self.pp,'fft_correlate_images','fft_ms'),(self.pp,'correlation_to_displacement','peak_ms'),(self.pp,'sig2noise_ratio','snr_ms')]
  elif self.name=='B3':targets=[(self.tp,'correalte_fft','fft_ms'),(self.tp,'correlation_to_displacement','peak_and_validation_and_download_ms')]
  elif self.name=='B4':targets=[(self.gp.CorrelationGPU,'correlate_windows','fft_ms'),(self.gp.CorrelationGPU,'get_subpixel_peak','peak_ms'),(self.gp.PIVGPU,'validate_fields','native_validation_ms')]
  try:
   for obj,name,key in targets:
    old=getattr(obj,name);original.append((obj,name,old))
    def wrapped(*a,_old=old,_key=key,**kw):
     self.sync();t=time.perf_counter();result=_old(*a,**kw);self.sync();self.parts[_key]=self.parts.get(_key,0)+(time.perf_counter()-t)*1000;return result
    setattr(obj,name,wrapped)
   yield
  finally:
   for obj,name,old in original:setattr(obj,name,old)
 def core(self,a,b,selection=None,actual_coordinates=False):
  ws,ov=self.ws,self.ov;native=None;snr=None;upload=download=None
  if selection is not None and self.name in ['B0','B1','B2','B3']:
   start=time.perf_counter();key=('packed',ws)
   if key not in self.cache:self.cache[key]=Backend(self.name,ws,0,profile=self.profile)
   # Gather official windows into a non-overlapping strip. The unchanged public
   # solver then treats each strip cell as one independent interrogation window.
   kw={'backend':'rust' if self.name=='B2' else 'scipy'} if 'backend' in inspect.signature(pyprocess.sliding_window_array).parameters else {}
   aa=pyprocess.sliding_window_array(a,ws,ov,**kw)[selection.ravel()].reshape(-1,ws)
   bb=pyprocess.sliding_window_array(b,ws,ov,**kw)[selection.ravel()].reshape(-1,ws)
   p=self.cache[key].core(aa,bb);x,y=pyprocess.get_coordinates(a.shape,ws,ov,center_on_field=False)
   mapped={}
   for k in ['u','v','snr']:
    mapped[k]=np.full(selection.shape,np.nan)
    if p[k] is not None:mapped[k][selection]=p[k].ravel()
    elif k=='snr':mapped[k]=None
   native=np.ones(selection.shape,bool);native[selection]=p['native'].ravel() if p['native'] is not None else False
   return dict(x=x,y=y,**mapped,native=native,core_ms=(time.perf_counter()-start)*1000,components_ms=p['components_ms'],upload_ms=p['upload_ms'],download_ms=p['download_ms'])
  self.sync();start=time.perf_counter()
  with self.instrument(),warnings.catch_warnings():
   warnings.simplefilter('ignore',RuntimeWarning)
   if self.name in ['B0','B1','B2']:
    kw={} if self.name=='B0' else {'backend':'scipy' if self.name=='B1' else 'rust'}
    u,v,snr=self.pp.extended_search_area_piv(a,b,ws,overlap=ov,dt=1,subpixel_method='gaussian',sig2noise_method='peak2peak',normalized_correlation=True,**kw)
    x,y=self.pp.get_coordinates(a.shape,ws,ov)
   elif self.name=='B3':
    t=time.perf_counter();aa=self.torch.as_tensor(a,device='cuda');bb=self.torch.as_tensor(b,device='cuda');self.sync();upload=(time.perf_counter()-t)*1000
    u,v,x,y,native=self.tp.extended_search_area_piv(aa,bb,window_size=ws,overlap=ov,validate=True,validation_ratio=1.2)
   elif self.name=='B4':
    key=(a.shape,ws,ov)
    if key not in self.cache:self.cache[key]=self.gp.piv_gpu(a.shape,ws,overlap_ratio=float(ov/ws),search_size_iters=1,n_fft=1,s2n_tol=1.2,median_tol=None,mad_tol=None,mean_tol=None,rms_tol=None,num_replacing_iters=0,smooth=False,center=False,dtype_f='float32')
    solver=self.cache[key]
    solver.gpu_process.piv_fields[0].field_mask=self.cp.asarray(~selection) if selection is not None else None
    u,v=solver(a,b);x,y=solver.coords
    # Upstream s2n_ratio property returns None on a second access; use its
    # already-computed public correlation attribute, without changing kernels.
    snr=solver.gpu_process.corr.sig2noise.get();native=solver.val_locations|solver.outliers
   else:raise ValueError(self.name)
  if actual_coordinates and self.name in ['B0','B1','B2','B3']:
   # Upstream window stacks start at pixel zero; get_coordinates historically
   # adds a centering offset even when no matching image crop was applied.
   x,y=pyprocess.get_coordinates(a.shape,ws,ov,center_on_field=False)
  self.sync();core_ms=(time.perf_counter()-start)*1000
  return dict(x=x,y=y,u=u,v=v,snr=snr,native=native,core_ms=core_ms,components_ms=dict(self.parts),upload_ms=upload,download_ms=download)
 def run(self,a,b,dt,mode='local',background_a=None,background_b=None,do_fit=True,selection=None,actual_coordinates=False):
  start=time.perf_counter();pa=preprocess(a,mode,background_a);pb=preprocess(b,mode,background_b);pre_ms=(time.perf_counter()-start)*1000
  r=self.core(pa,pb,selection=selection,actual_coordinates=actual_coordinates);x,y,u,v=r['x'],r['y'],r['u'],r['v'];t=time.perf_counter()
  bad=~np.isfinite(u+v)
  if r['snr'] is not None:bad|=~np.isfinite(r['snr'])|(r['snr']<1.2)
  if r['native'] is not None:bad|=r['native']
  if selection is not None:bad|=~selection
  kw={'backend':'rust' if self.name=='B2' else 'scipy'} if 'backend' in inspect.signature(validation.local_median_val).parameters else {}
  spatial_mask=~selection if selection is not None else np.zeros_like(bad)
  bad|=validation.local_median_val(np.ma.array(u,mask=spatial_mask),np.ma.array(v,mask=spatial_mask),2.5,2.5,size=1,**kw)
  texture=np.minimum(window_std(pa,x,y,self.ws),window_std(pb,x,y,self.ws));bad|=texture<2.
  validation_ms=(time.perf_counter()-t)*1000;t=time.perf_counter();good=~bad
  est=fit_affine(np.c_[x[good],y[good]],np.c_[u[good],v[good]],dt,a.shape,self.name) if do_fit else None
  fit_ms=(time.perf_counter()-t)*1000
  r.update(valid=good,texture=texture,estimate=est,preprocess_ms=pre_ms,validation_ms=validation_ms,fit_ms=fit_ms,total_ms=(time.perf_counter()-start)*1000,dtype=str(u.dtype),replacement_policy='NO_REPLACED_VECTORS_PUBLISHED',quality_scope='native TorchPIV ratio flag; ratio value unavailable' if self.name=='B3' else 'peak2peak')
  return r
