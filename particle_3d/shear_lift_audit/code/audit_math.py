"""Literature-checked, post-processing-only quantities. No trajectory integrator."""
import numpy as np

COEFFICIENT_RADIUS = 6.46
FORMULA_ROLE = 'ORDER_OF_MAGNITUDE_ONLY'

def local_tensors(gradient):
    g = np.asarray(gradient, float)
    e = (g + g.swapaxes(-1, -2)) / 2
    w = (g - g.swapaxes(-1, -2)) / 2
    curl = np.stack((g[...,2,1]-g[...,1,2],g[...,0,2]-g[...,2,0],g[...,1,0]-g[...,0,1]),axis=-1)
    return e,w,curl,np.sqrt(2*np.sum(e*e,axis=(-2,-1)))

def candidate(a, slip_speed, shear, mu, rho):
    a,s,g = np.broadcast_arrays(np.asarray(a,float),np.asarray(slip_speed,float),np.asarray(shear,float))
    if np.any(a<=0) or np.any(s<0) or np.any(g<0) or mu<=0 or rho<=0:
        raise ValueError('Positive size/fluid parameters and nonnegative speed/shear required')
    return COEFFICIENT_RADIUS*mu*a*a*s*np.sqrt(g/(mu/rho))

def dimensionless(a, slip, shear, gap, nu):
    a,s,g,h=np.broadcast_arrays(a,slip,shear,gap)
    rp=a*s/nu;rg=a*a*g/nu;root=np.sqrt(rg)
    margin=np.divide(rp,root,out=np.full(rp.shape,np.inf),where=root>0)
    lg=np.sqrt(np.divide(nu,g,out=np.full(g.shape,np.inf),where=g>0))
    return dict(Re_p=rp,Re_G=rg,h_over_a=h/a,slip_margin=margin,shear_margin=root,
                shear_length_m=lg,wall_margin=(a+h)/lg,oseen_length_m=np.divide(nu,s,out=np.full(s.shape,np.inf),where=s>0))

def validity(rp,rg,wall_margin, separation=.1):
    """Descriptive asymptotic screen, NOT a validated threshold or shell model."""
    root=np.sqrt(rg)
    return (root>0)&(root<separation)&(rp<separation*root)&(wall_margin>1/separation)

def force_ratio(numerator,denominator,zero_scale,numerator_zero_scale=None):
    n,d,t=np.broadcast_arrays(numerator,denominator,zero_scale)
    if np.any(n<0) or np.any(d<0) or np.any(t<0):raise ValueError('Magnitude only')
    ratio=np.full(n.shape,np.nan);good=d>t
    np.divide(n,d,out=ratio,where=good)
    status=np.full(n.shape,2,np.uint8) # DENOMINATOR_NEAR_ZERO_NOT_EVALUABLE
    nt=t if numerator_zero_scale is None else np.broadcast_to(numerator_zero_scale,n.shape)
    status[(n<=nt)&(d<=t)]=1 # BOTH_FORCE_MAGNITUDES_NEAR_ZERO
    status[good]=0
    return ratio,status

def wall_coefficient(a,h,mu):
    chi=np.asarray(h)/a;s=np.clip((chi-.01)/.04,0,1)
    weight=(1-s)**2*(1+2*s)
    return weight*6*np.pi*mu*a*a/np.maximum(h,np.maximum(2e-9,1e-3*a))

def align_saved(samples):
    """Saved endpoint velocity belongs to previous accepted point; keep endpoint too."""
    s=np.asarray(samples)
    x=s[:,1:4].copy();x[1:]=s[:-1,1:4]
    if len(s)>1:
        if not np.allclose(s[1:,19],s[:-1,0],rtol=0,atol=2e-14):
            raise ValueError('Unsupported saved velocity evaluation clock')
        if not np.allclose(s[1:,0]-s[:-1,0],s[1:,18],rtol=1e-10,atol=2e-14):
            raise ValueError('Saved step duration mismatch')
        if not np.allclose(s[1:,1:4],x[1:]+s[1:,18,None]*s[1:,4:7],rtol=0,atol=5e-19):
            raise ValueError('Saved position/held-velocity identity mismatch')
    return x

def stats(values,weights=None):
    a=np.asarray(values);good=np.isfinite(a)
    if weights is not None:good&=np.asarray(weights)>0
    x=a[good]
    if not len(x):return dict(count=0,median=None,P90=None,P95=None,P99=None,max=None)
    if weights is None:q=np.quantile(x,[.5,.9,.95,.99])
    else:
        order=np.argsort(x);x=x[order];w=np.asarray(weights)[good][order]
        q=np.interp([.5,.9,.95,.99],np.cumsum(w)/np.sum(w),x)
    return dict(count=int(len(x)),**dict(zip(['median','P90','P95','P99'],map(float,q))),max=float(np.max(x)))
