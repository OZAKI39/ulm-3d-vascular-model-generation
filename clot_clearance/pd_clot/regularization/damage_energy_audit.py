"""Read-only stored-energy diagnostics. NOT a damage law or mechanics solver.

All force/energy evaluations reuse the existing production kernel. Differences
are signed old-minus-new. No clipping, absolute-value repair or time advance.
"""
from dataclasses import dataclass
import numpy as np
from ..mechanics import shape_data
from ..fragment_mechanics import prepare_supported_shape
from .fast_mechanics import evaluate_supported
from .energy import potential
from ..fragment_topology import surface_particles


@dataclass
class EnergyState:
    particle: dict
    global_energy: dict
    shape: tuple
    force_N: np.ndarray


def audit_state(cloud, x, material, safety, *, g=None, D=None, frozen_shape=None):
    if (g is None) == (D is None): raise ValueError('Supply exactly one of g or D')
    integrity=np.array(g if g is not None else 1-np.asarray(D),float,copy=True)
    coords=np.array(x,float,copy=True)
    if coords.shape!=cloud.X.shape or not np.isfinite(coords).all():raise ValueError('Invalid coordinates')
    if integrity.shape!=(len(cloud.pairs),) or not np.isfinite(integrity).all() or np.any((integrity<0)|(integrity>1)):
        raise ValueError('Integrity must be finite in [0,1]')
    shape=prepare_supported_shape(cloud,integrity,safety)
    K,degree,remaining=shape_data(len(cloud.X),cloud.pairs,cloud.xi,cloud.weight,integrity,cloud.volume)
    eig=np.linalg.eigvalsh(K)
    condition=np.full(len(cloud.X),np.nan)
    for i,rank in enumerate(shape[3]):
        if rank:condition[i]=eig[i,-1]/eig[i,3-rank]
    evaluation_shape=shape
    if frozen_shape is not None:
        if not np.all(shape[3]==3) or not np.all(frozen_shape[3]==3):raise ValueError('Frozen-operator counterfactual is full-rank only')
        evaluation_shape=(frozen_shape[0],frozen_shape[1],remaining,shape[3],frozen_shape[4])
    m=evaluate_supported(cloud,coords,integrity,material,safety,evaluation_shape)
    constitutive,stabilization=potential(cloud,m)
    W=np.divide(m[3],remaining,out=np.zeros_like(remaining),where=remaining>0)
    p=dict(remaining=remaining,support_rank=shape[3],K=K,K_eigenvalues=eig,K_condition=condition,
           degree=degree,F=m[1],J=m[2],W_undegraded_Pa=W,constitutive_energy_density=m[3],
           stabilization_energy_density=m[4],constitutive_energy_J=cloud.volume*m[3],
           stabilization_energy_J=cloud.volume*m[4],total_energy_J=cloud.volume*(m[3]+m[4]))
    return EnergyState(p,dict(U_constitutive_J=constitutive,U_stabilization_J=stabilization,
                             U_total_J=constitutive+stabilization),shape,m[0])


def decompose(cloud, old, new):
    a,b=old.particle,new.particle;V=cloud.volume
    soft=V*(a['remaining']-b['remaining'])*a['W_undegraded_Pa']
    reconstruction=V*b['remaining']*(a['W_undegraded_Pa']-b['W_undegraded_Pa'])
    constitutive=a['constitutive_energy_J']-b['constitutive_energy_J']
    stabilization=a['stabilization_energy_J']-b['stabilization_energy_J']
    total=constitutive+stabilization
    norm=lambda A:np.linalg.norm(A,axis=(1,2))
    dK=norm(b['K']-a['K']);dF=norm(b['F']-a['F'])
    p=dict(delta_U_softening=soft,delta_U_reconstruction=reconstruction,
           delta_U_stabilization=stabilization,delta_U_constitutive=constitutive,delta_U_total=total,
           decomposition_residual=total-soft-reconstruction-stabilization,
           absolute_K_change=dK,relative_K_change=dK/np.maximum(norm(a['K']),np.finfo(float).tiny),
           absolute_F_change=dF,relative_F_change=dF/np.maximum(norm(a['F']),np.finfo(float).tiny),
           K_condition_change=b['K_condition']-a['K_condition'],
           K_anisotropy_change=(b['K_eigenvalues'][:,-1]-b['K_eigenvalues'][:,0])/np.maximum(b['K_eigenvalues'].sum(axis=1),np.finfo(float).tiny)-
                               (a['K_eigenvalues'][:,-1]-a['K_eigenvalues'][:,0])/np.maximum(a['K_eigenvalues'].sum(axis=1),np.finfo(float).tiny))
    direct=old.global_energy['U_total_J']-new.global_energy['U_total_J']
    s={key+'_J':float(value.sum()) for key,value in p.items() if key.startswith('delta_U_')}
    s['delta_U_total_J']=direct
    s.update(decomposition_residual_J=direct-s['delta_U_softening_J']-s['delta_U_reconstruction_J']-s['delta_U_stabilization_J'],
             maximum_particle_decomposition_residual_J=float(np.max(abs(p['decomposition_residual']))),
             max_relative_F_change=float(p['relative_F_change'].max()),max_relative_K_change=float(p['relative_K_change'].max()),
             max_K_condition=float(np.nanmax(b['K_condition'])))
    return p,s


def thermodynamic_force(cloud,x,D,bond_ids,material,safety,epsilons,relative_tolerance=.03,absolute_tolerance=1e-22):
    """Y=-dU/dD at fixed x; full epsilon history, explicit unresolved signs.

    Difference per-particle energies BEFORE summing. This avoids cancellation
    against energies of unaffected particles without altering kernel values.
    """
    D=np.array(D,float,copy=True);epsilons=np.array(epsilons,float)
    if np.any(epsilons<=0) or np.any(np.diff(epsilons)<=0):raise ValueError('Increasing positive epsilon list required')
    if np.any((D<0)|(D>1)):raise ValueError('Invalid D')
    ids=np.array(bond_ids,int);baseline=audit_state(cloud,x,material,safety,D=D)
    Y=np.empty((len(ids),len(epsilons)));soft=np.empty_like(Y);rec=np.empty_like(Y);stab=np.empty_like(Y);dF=np.empty_like(Y)
    schemes=np.empty(Y.shape,dtype='U8')
    for i,k in enumerate(ids):
        for j,e in enumerate(epsilons):
            low=D.copy();high=D.copy()
            if D[k]>=e and D[k]+e<=1:
                low[k]-=e;high[k]+=e;denom=2*e;scheme='central'
                a=audit_state(cloud,x,material,safety,D=low);b=audit_state(cloud,x,material,safety,D=high)
            elif D[k]+e<=1:
                high[k]+=e;denom=e;scheme='forward';a=baseline;b=audit_state(cloud,x,material,safety,D=high)
            elif D[k]>=e:
                low[k]-=e;denom=e;scheme='backward';a=audit_state(cloud,x,material,safety,D=low);b=baseline
            else:raise ValueError('Epsilon too large for bounded damage')
            p,_=decompose(cloud,a,b)
            Y[i,j]=p['delta_U_total'].sum()/denom
            soft[i,j]=p['delta_U_softening'].sum()/denom;rec[i,j]=p['delta_U_reconstruction'].sum()/denom
            stab[i,j]=p['delta_U_stabilization'].sum()/denom
            dF[i,j]=float(np.linalg.norm(b.particle['F']-a.particle['F'])/denom);schemes[i,j]=scheme
    estimate=np.full(len(ids),np.nan);error=np.full(len(ids),np.nan);plateau=np.full(len(ids),-1,int)
    stable=np.zeros(len(ids),bool);resolved=np.zeros(len(ids),bool)
    for i in range(len(ids)):
        for j in range(len(epsilons)-2):
            samples=Y[i,j:j+3];median=float(np.median(samples));spread=float(np.ptp(samples))
            if spread<=absolute_tolerance+relative_tolerance*abs(median):
                estimate[i]=median;error[i]=spread+absolute_tolerance;plateau[i]=j;stable[i]=True
                resolved[i]=abs(median)>error[i];break
    return dict(bond_ids=ids,epsilons=epsilons,Y_J=Y,Y_softening_J=soft,Y_reconstruction_J=rec,Y_stabilization_J=stab,
                F_change_per_damage=dF,scheme=schemes,estimate_J=estimate,uncertainty_J=error,plateau_start=plateau,
                stable=stable,sign_resolved=resolved,negative_resolved=resolved & (estimate<0))


def categories(cloud,c,g):
    surface,support=surface_particles(cloud,g,c['surface'])
    fixed_near=cloud.X[:,2]<=cloud.X[cloud.fixed,2].max()+cloud.spacing*(1+1e-10)
    distance=np.linalg.norm(cloud.X-np.asarray(c['streaming']['bubble_center_m']),axis=1)
    source_near=distance<=c['localization']['influence_radius_m']
    interior=~surface & ~fixed_near & ~source_near
    masks=dict(INTERIOR=interior,FREE_SURFACE=surface,FIXED_BASE_NEAR=fixed_near,SOURCE_NEAR=source_near)
    label=np.full(len(cloud.X),'OTHER',dtype='U20')
    for name in ['INTERIOR','FREE_SURFACE','SOURCE_NEAR','FIXED_BASE_NEAR']:label[masks[name]]=name
    return dict(**masks,OTHER=label=='OTHER',exclusive_category=label,support_fraction=support,
                distance_to_source_m=distance,distance_to_fixed_layer_m=np.maximum(cloud.X[:,2]-cloud.X[cloud.fixed,2].max(),0.))


def counterfactuals(cloud,x,g0,g1,material,safety):
    a=audit_state(cloud,x,material,safety,g=g0);b=audit_state(cloud,x,material,safety,g=g1)
    frozen_const=float(np.sum(cloud.volume*b.particle['remaining']*a.particle['W_undegraded_Pa']))
    frozen_operator=audit_state(cloud,x,material,safety,g=g1,frozen_shape=a.shape)
    return dict(label='DIAGNOSTIC ONLY / NOT A CANDIDATE MODEL / NOT USED IN PRODUCTION',
        actual_delta_U_J=a.global_energy['U_total_J']-b.global_energy['U_total_J'],
        freeze_F_constitutive_only_delta_U_J=a.global_energy['U_constitutive_J']-frozen_const,
        freeze_old_inverse_and_trace_delta_U_J=a.global_energy['U_total_J']-frozen_operator.global_energy['U_total_J'],
        constitutive_only_delta_U_J=a.global_energy['U_constitutive_J']-b.global_energy['U_constitutive_J'],
        stabilization_only_delta_U_J=a.global_energy['U_stabilization_J']-b.global_energy['U_stabilization_J'],
        frozen_K_meaning='New damaged numerator A and remaining; old K inverse and trace; not affine reproducing',
        freeze_F_meaning='Old F/W; new remaining in constitutive term; old stabilization held unchanged')
