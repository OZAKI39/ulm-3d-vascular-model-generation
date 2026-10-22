import time, warnings
import numpy as np
import cv2
from openpiv import pyprocess, validation, windef, filters
from pycpd import AffineRegistration
from common import P,NAMES
from affine import fit_affine,parameters,empty

def preprocess(gray, mode, objects=None, background=None):
    gray=gray.astype(np.uint8)
    if mode=='raw':return gray.astype(np.float32)
    signal=cv2.subtract(cv2.GaussianBlur(gray,(0,0),4),gray).astype(np.float32)
    if mode=='local':return signal
    if mode=='fixed_suppressed':
        # Causal temporal background, built from frames available BEFORE endpoint processing.
        return np.maximum(np.asarray(background,np.float32)-gray,0) if background is not None else signal
    if mode=='candidate':
        m=np.zeros_like(gray)
        for o in objects or []:
            x,y=o['xy'];r=max(4,min(16,int(np.sqrt(o['area']/np.pi)*2)))
            cv2.circle(m,(round(x),round(y)),r,1,-1)
        return signal*m
    raise ValueError(mode)

def piv(a,b,dt,mode='local',objects_a=None,objects_b=None,background_a=None,background_b=None,multipass=False,bootstrap=0):
    t=time.perf_counter();a=preprocess(a,mode,objects_a,background_a);b=preprocess(b,mode,objects_b,background_b)
    ws=P['piv']['window_size'];ov=P['piv']['overlap']
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        u,v,snr=pyprocess.extended_search_area_piv(a,b,ws,overlap=ov,dt=1,subpixel_method='gaussian',sig2noise_method='peak2peak',normalized_correlation=True)
    x,y=pyprocess.get_coordinates(a.shape,ws,ov)
    bad=~np.isfinite(u+v+snr)|(snr<P['piv']['snr_min'])
    bad|=validation.local_median_val(np.ma.array(u,mask=np.zeros_like(bad)),np.ma.array(v,mask=np.zeros_like(bad)),2.5,2.5,size=1)
    coarse_peak_quality=None
    if multipass:
        # Official multipass implementation; retain final flags and independently
        # reject weak original-image texture. Never publish its replaced outliers.
        settings=windef.PIVSettings();settings.windowsizes=(64,32);settings.overlap=(32,16);settings.num_iterations=2
        settings.sig2noise_method='peak2peak';settings.sig2noise_threshold=P['piv']['snr_min'];settings.sig2noise_validate=True
        settings.show_all_plots=False
        try:
            x0,y0,u0,v0,s0=windef.first_pass(a,b,settings)
            coarse_peak_quality=float(np.nanmedian(s0))
            u0=np.ma.array(u0,mask=np.zeros_like(u0,dtype=bool))
            v0=np.ma.array(v0,mask=np.zeros_like(v0,dtype=bool))
            flags=validation.typical_validation(u0,v0,s0,settings)
            if np.all(flags):raise ValueError('all coarse windows invalid')
            u0,v0=filters.replace_outliers(np.ma.array(u0),np.ma.array(v0),flags)
            x,y,u,v,grid,bad=windef.multipass_img_deform(a,b,1,x0,y0,u0,v0,settings)
            u=np.asarray(u);v=np.asarray(v);bad=np.asarray(bad)|grid|~np.isfinite(u+v)
            snr=None;ws=32
        except Exception as e:
            return {'estimate':empty('MULTIPASS_FAILED:'+str(e),'OpenPIV_multipass'),'algorithm_ms':(time.perf_counter()-t)*1000,'vectors':[],'valid_vector_count':0,'rejected_fraction':1.}
    texture=[]
    for xx,yy in zip(x.ravel(),y.ravel()):
        xi,yi=int(xx),int(yy);r=ws//2
        texture.append(min(np.std(a[yi-r:yi+r,xi-r:xi+r]),np.std(b[yi-r:yi+r,xi-r:xi+r])))
    texture=np.array(texture).reshape(x.shape);bad|=texture<P['piv']['texture_std_min']
    good=~bad;algorithm_ms=(time.perf_counter()-t)*1000
    est=fit_affine(np.c_[x[good],y[good]],np.c_[u[good],v[good]],dt,a.shape,'OpenPIV_'+('multipass_' if multipass else '')+mode,bootstrap)
    vec=np.c_[x.ravel(),y.ravel(),u.ravel(),v.ravel(),good.ravel().astype(float)]
    return {'estimate':est,'algorithm_ms':algorithm_ms,'vectors':vec.tolist(),'valid_vector_count':int(good.sum()),'rejected_fraction':float(bad.mean()),'peak_quality_median':float(np.nanmedian(snr)) if snr is not None else None,'coarse_peak_quality_median':coarse_peak_quality,'peak_quality_scope':'final flags only for multipass; no fabricated final SNR','texture_std_median':float(np.median(texture))}

def cpd(x,y,dt,shape,multistart=True):
    t0=time.perf_counter();x=np.asarray(x,float).reshape(-1,2);y=np.asarray(y,float).reshape(-1,2)
    est=empty('INSUFFICIENT_POINTS','pycpd_affine')
    if min(len(x),len(y))<P['cpd']['min_points']:return {'estimate':est,'algorithm_ms':(time.perf_counter()-t0)*1000}
    center=np.array([shape[1]/2,shape[0]/2]);scale=max(shape)/2
    xx=(x-center)/scale;yy=(y-center)/scale; fits=[]
    try:
        for theta in ([0,-.18,.18] if multistart else [0]):
            rot=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
            r=AffineRegistration(X=yy,Y=xx,B=rot,max_iterations=P['cpd']['max_iterations'],tolerance=P['cpd']['tolerance'],w=P['cpd']['outlier_weight'])
            transformed,(B,tt)=r.register();F=B.T;b=center-F@center+tt.ravel()*scale
            p=parameters(F,b,dt,center)
            # Responsibility-weighted residual is available, but sigma2 is NOT truth accuracy.
            fits.append({'F':F,'b':b,'params':p,'sigma2':r.sigma2,'iteration':r.iteration,'diff':r.diff})
        best=min(fits,key=lambda z:z['sigma2']);ps=np.array([z['params'] for z in fits]);value=best['params']
        eig=np.linalg.eigvalsh(np.cov(x.T));ratio=np.sqrt(eig[-1]/max(eig[0],1e-12));change=abs(len(x)-len(y))/max(len(x),len(y))
        est.update(sample_count=len(x),F=best['F'].tolist(),b=best['b'].tolist(),raw_estimate=value.tolist(),residual=float(np.sqrt(best['sigma2'])*scale),iterations=best['iteration'],converged=best['diff']<P['cpd']['tolerance'],count_change_fraction=change,axis_ratio=ratio,
                   uncertainty={'method':'three_initial_rotation_empirical_range_not_CI','calibrated':False,'min':ps.min(0).tolist(),'max':ps.max(0).tolist()})
        reasons=[None]*6
        if change>P['cpd']['count_change_max']:reasons=['DETECTION_COUNT_CHANGE']*6
        elif not est['converged']:reasons=['CPD_ITERATION_LIMIT']*6
        elif est['residual']>3:reasons=['POOR_POINT_REGISTRATION']*6
        elif ratio<P['cpd']['isotropy_ratio_min']:reasons=['SYMMETRIC_POINTSET_NO_UNIQUE_MATERIAL_MAP']*6
        elif np.linalg.cond(np.c_[np.ones(len(xx)),xx])>P['piv']['condition_max']:reasons=['DEGENERATE_SPATIAL_SUPPORT']*6
        else:
            for j,z in enumerate(np.ptp(ps,axis=0)):
                if z>P['parameter_error_limits'][j]:reasons[j]='MULTISTART_AMBIGUITY'
        est['values']={k:float(value[i]) if reasons[i] is None else None for i,k in enumerate(NAMES)}
        est['valid']={k:reasons[i] is None for i,k in enumerate(NAMES)};est['reasons']=dict(zip(NAMES,reasons))
        est['observability_status']='INTERVAL_AFFINE_ESTIMATE' if all(v is None for v in reasons) else ('PARTIAL' if any(v is None for v in reasons) else 'NOT_OBSERVABLE')
    except Exception as e:est['failure']=repr(e)
    return {'estimate':est,'algorithm_ms':(time.perf_counter()-t0)*1000}
