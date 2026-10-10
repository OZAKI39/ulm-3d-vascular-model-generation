"""Fixed spatial field, sampled at current x; no adaptive force renormalization."""
import numpy as np

def spatial_shape(X,r):
    R=r['radius_m'];w=.4*R;rad=np.sqrt(np.sum(X*X,axis=1)+1e-30);n=X/rad[:,None]
    if r['pattern_id']=='LOCAL_GAUSSIAN_SHEAR':
        v=np.zeros_like(X);v[:,0]=np.exp(-np.sum((X-[0,0,R])**2,axis=1)/(2*w*w))
    elif r['pattern_id']=='OPPOSING_SHEAR':
        a=np.exp(-np.sum((X-[.5*R,0,.85*R])**2,axis=1)/(2*w*w))
        b=np.exp(-np.sum((X-[-.5*R,0,.85*R])**2,axis=1)/(2*w*w));v=np.zeros_like(X);v[:,0]=a-b
    elif r['pattern_id']=='VORTEX_TANGENTIAL':
        env=np.exp(-np.sum((X-[0,0,.5*R])**2,axis=1)/(2*R*R));v=np.c_[-X[:,1],X[:,0],np.zeros(len(X))]/R*env[:,None]
    else:raise ValueError(r['pattern_id'])
    return v-n*np.sum(v*n,axis=1)[:,None]

def pulse(t,T):
    q=np.asarray(t)/T
    return np.where(q<0,0,np.where(q<.25,.5*(1-np.cos(4*np.pi*q)),np.where(q<.75,1.,np.where(q<1,.5*(1+np.cos(4*np.pi*(q-.75))),0.))))

def prepare_provider(g,r):
    shape=spatial_shape(g['X'],r);weighted=np.linalg.norm(shape,axis=1)*g['areas'];area=g['areas'].sum()
    factor=r['tau_budget_Pa']*area/weighted.sum();F=shape*g['areas'][:,None]*factor
    return dict(scale_Pa=factor,reference_force_N=F,reference_peak_traction_Pa=float(np.linalg.norm(shape[g['surface']],axis=1).max()*factor),
        initial_integrated_force_budget_N=float(np.linalg.norm(F,axis=1).sum()),reference_surface_area_m2=float(area),
        normalization='Equal reference-geometry plateau mean traction; frozen scale. Actual budgets change as bodies move; no renormalization that defeats spatial decay',
        wetted_set='Initial exposed curved surface only; new fracture surfaces not added',direction='Tangent to fixed Eulerian radial direction at current position, not a tracked deformed surface normal')

def lammps_provider(g,r,load):
    R=r['radius_m'];w=.4*R;T=r['duration_s'];f=load['scale_Pa'];lines=[]
    def atom(k,v):lines.append(f'variable {k} atom {v}')
    atom('rr','sqrt(x*x+y*y+z*z+1e-30)')
    for a in 'xyz':atom('n'+a,f'{a}/v_rr')
    if r['pattern_id']=='LOCAL_GAUSSIAN_SHEAR':atom('rawx',f'exp(-(x*x+y*y+(z-{R:.17g})^2)/{2*w*w:.17g})');atom('rawy','0.0')
    elif r['pattern_id']=='OPPOSING_SHEAR':
        atom('rawx',f'exp(-((x-{.5*R:.17g})^2+y*y+(z-{.85*R:.17g})^2)/{2*w*w:.17g})-exp(-((x+{.5*R:.17g})^2+y*y+(z-{.85*R:.17g})^2)/{2*w*w:.17g})');atom('rawy','0.0')
    else:
        atom('env',f'exp(-(x*x+y*y+(z-{.5*R:.17g})^2)/{2*R*R:.17g})');atom('rawx',f'-y/{R:.17g}*v_env');atom('rawy',f'x/{R:.17g}*v_env')
    atom('rawz','0.0');atom('dotn','v_rawx*v_nx+v_rawy*v_ny+v_rawz*v_nz')
    lines.append(f'variable amp equal ternary(time<{T/4:.17g},0.5*(1-cos(4*PI*time/{T:.17g})),ternary(time<{.75*T:.17g},1.0,ternary(time<{T:.17g},0.5*(1+cos(4*PI*(time/{T:.17g}-0.75))),0.0)))')
    for a in 'xyz':atom('load_f'+a,f'd_area*{f:.17g}*v_amp*(v_raw{a}-v_dotn*v_n{a})')
    lines.append('fix stream surface addforce v_load_fx v_load_fy v_load_fz')
    return '\n'.join(lines)+'\n'
