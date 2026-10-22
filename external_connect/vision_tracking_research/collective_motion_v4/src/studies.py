"""Known-motion scoring and official multipass calls; replacements stay flagged."""
import time,warnings
import numpy as np
from backends import preprocess,window_std
from xcommon import fit_affine,parameters,NAMES
from openpiv import windef,validation,filters
def multipass(engine,a,b,dt,mode='local',background_a=None,background_b=None):
 start=time.perf_counter();aa=preprocess(a,mode,background_a);bb=preprocess(b,mode,background_b);ws=32
 if engine.name=='B4':
  key=('multi',a.shape)
  if key not in engine.cache:engine.cache[key]=engine.gp.piv_gpu(a.shape,32,search_size_iters=(1,1),overlap_ratio=.5,n_fft=2,s2n_tol=1.2,median_tol=None,num_replacing_iters=(2,0),replacing_method='median',smooth=False,center=False)
  sol=engine.cache[key];u,v=sol(aa,bb);x,y=sol.coords;snr=sol.gpu_process.corr.sig2noise.get();bad=sol.val_locations|sol.outliers
 elif engine.name in ['B1','B2']:
  settings=windef.PIVSettings();settings.windowsizes=(64,32);settings.overlap=(32,16);settings.num_iterations=2;settings.backend='rust' if engine.name=='B2' else 'scipy';settings.sig2noise_method='peak2peak';settings.sig2noise_threshold=1.2;settings.sig2noise_validate=True;settings.show_all_plots=False
  x0,y0,u0,v0,s0=windef.first_pass(aa,bb,settings);u0=np.ma.array(u0,mask=np.zeros_like(u0,dtype=bool));v0=np.ma.array(v0,mask=np.zeros_like(v0,dtype=bool));flags=validation.typical_validation(u0,v0,s0,settings)
  if np.all(flags):raise ValueError('ALL_COARSE_VECTORS_INVALID')
  u0,v0=filters.replace_outliers(u0,v0,flags)
  x,y,u,v,grid,bad=windef.multipass_img_deform(aa,bb,1,x0,y0,u0,v0,settings);u=np.asarray(u);v=np.asarray(v);bad=np.asarray(bad)|grid;snr=None
 elif engine.name=='B3':
  ta=engine.torch.as_tensor(aa,device='cuda');tb=engine.torch.as_tensor(bb,device='cuda')
  u,v,x,y,bad=engine.tp.extended_search_area_piv(ta,tb,64,32,True,1.2)
  key=('multi',a.shape)
  if key not in engine.cache:engine.cache[key]=engine.tp.piv_iteration_CWS(a.shape,32,16,engine.torch.device('cuda'))
  u,v,x,y,bad=engine.cache[key](ta,tb,x,y,u,v,bad);snr=None
 else:raise ValueError(engine.name)
 engine.sync();core_ms=(time.perf_counter()-start)*1000
 bad=np.asarray(bad)|~np.isfinite(u+v)
 if snr is not None:bad|=~np.isfinite(snr)|(snr<1.2)
 bad|=validation.local_median_val(np.ma.array(u),np.ma.array(v),2.5,2.5,size=1)
 bad|=np.minimum(window_std(aa,x,y,ws),window_std(bb,x,y,ws))<2
 est=fit_affine(np.c_[x[~bad],y[~bad]],np.c_[u[~bad],v[~bad]],dt,a.shape,engine.name+'_MULTIPASS')
 return dict(x=x,y=y,u=u,v=v,valid=~bad,snr=snr,estimate=est,total_ms=(time.perf_counter()-start)*1000,core_ms=core_ms,interpolated_predictor=True,replacement_policy='COARSE_PREDICTOR_REPLACEMENT_ONLY_FINAL_NATIVE_FLAGS_REJECTED',final_window=32)
def score(r,data,meta,start,end):
 dt=float(data['times'][end]-data['times'][start]);H=data['H'][end]@np.linalg.inv(data['H'][start]);xy=np.c_[r['x'].ravel(),r['y'].ravel()];expected=(np.c_[xy,np.ones(len(xy))]@H.T)[:,:2]-xy
 if meta['case']=='local_opposite':expected[:,0]+=np.where(xy[:,1]<128,26.,-26.)*dt
 observed=np.c_[r['u'].ravel(),r['v'].ravel()];epe=np.linalg.norm(observed-expected,axis=1);good=r['valid'].ravel();est=r['estimate'];accepted=all(est['valid'].values());truth=parameters(H[:2,:2],H[:2,2],dt,np.array([128.,128.])) if meta['known_global_affine'] else None
 err=np.abs(np.asarray(est['raw_estimate'])-truth) if truth is not None and est.get('raw_estimate') is not None else None
 badtruth=not meta['observable_truth'] or (meta['case']=='local_opposite')
 over=bool(np.any(err>np.array([5,5,.2,.2,.1,.1]))) if err is not None else None
 return dict(accepted=accepted,false_accept=bool(accepted and (badtruth or over)),correct_reject=bool(not accepted and badtruth),valid_vectors=int(good.sum()),n_vectors=len(good),epe_all=float(np.nanmean(epe)) if np.isfinite(epe).any() else None,epe_valid=float(np.mean(epe[good])) if good.any() else None,six_truth=truth,six_absolute_error=err,six_sign_correct=(np.sign(est['raw_estimate'])==np.sign(truth)).tolist() if err is not None else None,over_tolerance=over,estimate=est,total_ms=r['total_ms'],dt_s=dt,measurement_midpoint_age_ms=dt*500,nonaffine_truth=not meta['known_global_affine'],low_texture_truth=not meta['observable_truth'])
