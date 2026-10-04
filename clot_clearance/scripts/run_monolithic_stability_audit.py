"""Joint energy references, curvature and diagnostic substeps. Never a solver fix."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import copy
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigh, solve
from pd_clot.pfpd_candidate.coupling_stability_audit import (
    ROOT,OLD,OUT,AuditProblem,audit_config,physics_config,read_npz,read_json,save_json,
    load_audit_state,load_saved_iteration,deterministic_starts,green_minimum)
from pd_clot.pfpd_candidate.at1_evolution import at1_constrained_update,HighDamageGate
from pd_clot.pfpd_candidate.adr_adapter_at1 import equilibrate_tensile
from scripts.run_coupling_stability_audit import worker_budget


def gradient_point(key):
    if key!='failed_load_initial':return load_audit_state(key)
    p,data=load_saved_iteration('loading_unloading',10,0)
    return p,dict(x=data['x_start'],s=p.lower.copy())


def gradient_audit(key):
    a=audit_config();p,data=gradient_point(key);z=p.pack(data['x'],data['s']);f0,g=p.objective_gradient(z)
    rng=np.random.default_rng(a['gradient_seed']);rows=[]
    for index in range(a['gradient_directions']):
        direction=rng.normal(size=len(z));direction/=np.linalg.norm(direction)
        at_lower=z[p.nx:]-p.lower<=a['active_s_tolerance']
        direction[p.nx:][at_lower]=abs(direction[p.nx:][at_lower])
        scheme='forward_second_order' if at_lower.any() else 'central'
        analytic=float(g@direction)
        for eps in a['gradient_epsilon']:
            fp=p.objective_gradient(z+eps*direction)[0]
            if scheme=='central':fd=(fp-p.objective_gradient(z-eps*direction)[0])/(2*eps)
            else:fd=(-3*f0+4*fp-p.objective_gradient(z+2*eps*direction)[0])/(2*eps)
            scale=max(abs(analytic),abs(fd),float(np.linalg.norm(g)),a['gradient_scale_floor'])
            probe_x,probe_s=p.unpack(z+eps*direction)
            rows.append(dict(direction=index,epsilon=eps,scheme=scheme,analytic=analytic,finite_difference=float(fd),
                normalized_error=float(abs(fd-analytic)/scale),scale=scale,
                probe_min_green=green_minimum(p.state(probe_x,probe_s))))
    passed=max(r['normalized_error'] for r in rows)<=a['gradient_relative_tolerance']
    save_json(OUT/f'monolithic/gradient_{key}.json',dict(passed=passed,metrics=p.metrics(data['x'],data['s']),
        energy_scale_J=p.energy_scale,displacement_scale_m=p.length_scale,rows=rows,
        explanation='Derivative of algebraic energy; out-of-domain probes are not physical candidates'))
    return dict(state=key,passed=passed,max_error=max(r['normalized_error'] for r in rows))


def full_hessian(problem,z,epsilon):
    _,g0=problem.objective_gradient(z);H=np.zeros((len(z),len(z)))
    schemes=[]
    for j in range(len(z)):
        unit=np.zeros_like(z);unit[j]=epsilon
        if j<problem.nx or (z[j]-problem.lower[j-problem.nx]>epsilon and problem.sc-z[j]>epsilon):
            H[:,j]=(problem.objective_gradient(z+unit)[1]-problem.objective_gradient(z-unit)[1])/(2*epsilon)
            schemes.append('central')
        elif problem.sc-z[j]>2*epsilon:
            H[:,j]=(-3*g0+4*problem.objective_gradient(z+unit)[1]-problem.objective_gradient(z+2*unit)[1])/(2*epsilon)
            schemes.append('forward_second_order')
        else:
            H[:,j]=(3*g0-4*problem.objective_gradient(z-unit)[1]+problem.objective_gradient(z-2*unit)[1])/(2*epsilon)
            schemes.append('backward_second_order')
    return H,schemes


def hessian_audit(task):
    label,x,s,strain,lower,epsilon=task
    a=audit_config();p=AuditProblem(physics_config(),strain,lower);z=p.pack(x,s)
    H,schemes=full_hessian(p,z,epsilon);scale=max(float(np.linalg.norm(H,2)),1.)
    symmetry=float(np.linalg.norm(H-H.T)/max(np.linalg.norm(H),1e-30));symmetric=.5*(H+H.T)
    state=p.state(x,s);expected=np.diag(2*state.W0*p.weights/p.energy_scale)
    hss_error=float(np.linalg.norm(H[p.nx:,p.nx:]-expected)/max(np.linalg.norm(expected),1e-30))
    groups=p.classify(x,s,a);interior=np.flatnonzero(groups['INTERIOR_ACTIVE'])
    keep=np.r_[np.arange(p.nx),p.nx+interior]
    reduced=symmetric[np.ix_(keep,keep)];evals,evecs=eigh(reduced)
    xx=symmetric[:p.nx,:p.nx];xs=symmetric[:p.nx,p.nx:][:,interior]
    ss=symmetric[p.nx:,p.nx:][np.ix_(interior,interior)]
    vals_x=eigh(xx,eigvals_only=True);cond=float(np.linalg.cond(xx))
    schur=None;schur_values=np.array([])
    if cond<a['hxx_condition_limit'] and vals_x[0]>a['curvature_relative_tolerance']*scale and len(interior):
        schur=ss-xs.T@solve(xx,xs,assume_a='pos');schur_values=eigh(.5*(schur+schur.T),eigvals_only=True)
    vector=np.zeros(len(z));vector[keep]=evecs[:,0]
    probes=[];base_value,gradient=p.objective_gradient(z)
    for amp in a['negative_mode_probe_amplitudes']:
        points=[]
        for sign in [-1,1]:
            trial=z+sign*amp*vector;tx,ts=p.unpack(trial)
            bound=bool(np.all(ts>=p.lower) and np.all(ts<p.sc))
            if not bound:points.append(dict(sign=sign,bounds=False));continue
            met=p.metrics(tx,ts);value=p.objective_gradient(trial)[0]
            points.append(dict(sign=sign,bounds=True,minimum_Green_strain=met['minimum_Green_strain'],
                in_tensile_domain=met['minimum_Green_strain']>=p.config['tolerances']['tensile_min_green_strain'],
                delta_normalized_energy=float(value-base_value),delta_energy_J=float((value-base_value)*p.energy_scale)))
        probes.append(dict(amplitude=amp,points=points))
    valid_matrix=symmetry<=a['hessian_symmetry_relative_tolerance'] and hss_error<=a['hss_relative_tolerance']
    threshold=a['curvature_relative_tolerance']*scale
    status=('NEGATIVE_CURVATURE_FOUND' if evals[0]<-threshold else
            'LOCALLY_POSITIVE' if evals[0]>threshold else 'INCONCLUSIVE') if valid_matrix else 'INCONCLUSIVE'
    folder=OUT/'hessian'/label/f'eps_{epsilon:g}';folder.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(folder/'BLOCKS.npz',H_raw=H,H_symmetric=symmetric,H_xx=H[:p.nx,:p.nx],
        H_xs=H[:p.nx,p.nx:],H_sx=H[p.nx:,:p.nx],H_ss=H[p.nx:,p.nx:],H_ss_analytic=expected,
        interior_phase_indices=interior,reduced=reduced,reduced_eigenvalues=evals,
        minimum_mode_scaled=vector,mechanical_eigenvalues=vals_x,
        Schur=np.empty((0,0)) if schur is None else schur,Schur_eigenvalues=schur_values,x=x,s=s,lower=lower)
    metrics=p.metrics(x,s)
    summary=dict(state=label,strain=strain,epsilon=epsilon,status=status,metrics=metrics,
        Hessian_symmetry_relative_error=symmetry,Hss_analytic_relative_error=hss_error,
        minimum_reduced_eigenvalue=float(evals[0]),minimum_mechanical_eigenvalue=float(vals_x[0]),
        minimum_Schur_eigenvalue=float(schur_values[0]) if len(schur_values) else None,
        Hxx_condition_number=cond,curvature_threshold=threshold,interior_phase_DOFs=len(interior),
        strong_lower_DOFs=int(groups['LOWER_BOUND_ACTIVE'].sum()),weak_lower_DOFs=int(groups['WEAK_LOWER_OR_UNCLASSIFIED'].sum()),
        scale=dict(displacement_m=p.length_scale,energy_J=p.energy_scale,phase=1.),
        gradient_dot_minimum_mode=float(gradient@vector),minimum_mode_probes=probes,
        no_symmetry_repair_before_error_measurement=True,schemes=schemes)
    save_json(folder/'SUMMARY.json',summary)
    return summary


class ReferenceExit(RuntimeError):
    def __init__(self,reason,z,metrics):super().__init__(reason);self.z=z.copy();self.metrics=metrics


def optimize_case(task):
    method,label,x,s,strain,lower=task;a=audit_config();p=AuditProblem(physics_config(),strain,lower)
    z0=p.pack(x,s);last_valid=z0.copy();evaluations=[];iterations=[]
    def fun(z):
        nonlocal last_valid
        xx,ss=p.unpack(z);met=p.metrics(xx,ss)
        if met['minimum_Green_strain']<p.config['tolerances']['tensile_min_green_strain']:
            raise ReferenceExit('MONOLITHIC_SOLUTION_LEFT_CURRENT_MODEL_DOMAIN',z,met)
        if met['maximum_s']>=p.sc:raise ReferenceExit('HIGH_DAMAGE_GATE_REACHED',z,met)
        value,gradient=p.objective_gradient(z);last_valid=z.copy()
        evaluations.append(dict(evaluation=len(evaluations),objective=value,**met))
        return value,gradient
    def callback(z):
        xx,ss=p.unpack(z);iterations.append(dict(iteration=len(iterations),**p.metrics(xx,ss)))
    bounds=[(None,None)]*p.nx+list(zip(p.lower,np.full(p.ns,p.sc)))
    exit_info=None;exception_info=None
    try:
        result=minimize(fun,z0,jac=True,method=method,bounds=bounds,callback=callback,options=a['optimizer_methods'][method])
        z=result.x;optimizer_success=bool(result.success);message=str(result.message);count=int(result.nit)
    except ReferenceExit as exc:
        z=last_valid.copy();optimizer_success=False;message=str(exc);count=len(iterations)
        exit_info=dict(reason=str(exc),metrics=exc.metrics);exit_z=exc.z.copy()
    except (ValueError,FloatingPointError) as exc:
        z=last_valid.copy();optimizer_success=False;message=str(exc);count=len(iterations)
        exception_info=dict(type=type(exc).__name__,message=str(exc))
    xx,ss=p.unpack(z);met=p.metrics(xx,ss)
    first_order=bool(optimizer_success and np.isfinite(met['total_J']) and
        met['minimum_lower_margin']>=-a['optimizer_bounds_tolerance'] and met['minimum_upper_margin']>0 and
        met['force_residual']<=p.config['mechanical_equilibrium']['relative_force_residual_tolerance'] and
        met['projected_phase_KKT']<=a['phase_KKT_tolerance'] and
        met['minimum_Green_strain']>=p.config['tolerances']['tensile_min_green_strain'])
    folder=OUT/'monolithic'/method/label;folder.mkdir(parents=True,exist_ok=True)
    arrays=dict(x=xx,s=ss,lower=lower,z=z,z_initial=z0)
    if exit_info:arrays['z_domain_or_sc_exit']=exit_z
    np.savez_compressed(folder/'STATE.npz',**arrays)
    summary=dict(method=method,start=label,strain=strain,optimizer_success=optimizer_success,
        message=message,iterations=count,evaluations=len(evaluations),first_order_candidate_passed=first_order,
        accepted_as_local_minimum=False,curvature_audit_pending=first_order,
        metrics=met,exit=exit_info,exception=exception_info,
        initial_metrics=p.metrics(x,s),maximum_x_distance_from_start_m=float(abs(xx-x).max()),
        maximum_s_distance_from_start=float(abs(ss-s).max()))
    save_json(folder/'SUMMARY.json',summary);save_json(folder/'HISTORY.json',dict(evaluations=evaluations,iterations=iterations))
    return summary


def subdivision_case(count):
    """Continue the saved load-009 state with unchanged existing subsolvers.

    This diagnostic driver duplicates no ADR/AT1 algorithm and enforces the
    original simultaneous convergence and domain gates verbatim.
    """
    a=audit_config();c=physics_config();prev=read_npz(OLD/'loading_unloading/load_009/accepted_state.npz')
    x=prev['x'].copy();s=prev['s_pf_bond'].copy();accepted=[];history=[];status='ENDPOINT_REACHED';reason=None
    folder=OUT/'continuation'/f'substeps_{count}';folder.mkdir(parents=True,exist_ok=False)
    loads=np.linspace(a['subdivision_start'],a['subdivision_endpoint'],count+1)[1:]
    save_json(folder/'CONFIG.json',dict(label='STEP_SUBDIVISION_DIAGNOSTIC_ONLY',loads=loads.tolist(),physics_config=c))
    for step,strain in enumerate(loads):
        p=AuditProblem(c,float(strain),s.copy());last_energy=None;converged=False
        for iteration in range(c['staggered']['maximum_iterations']):
            answer=p.phi(s,x)
            np.savez_compressed(folder/f'step_{step:02d}_iter_{iteration:03d}.npz',x=answer['x_equilibrium'],
                                s_before=s,s_after=answer.get('s_output',s),Y=answer['Y'])
            if answer['status']!='OK':status='STOPPED';reason=answer['status'];break
            x=answer['x_equilibrium'];new=answer['s_output'];met=answer['after']
            ds=float(abs(new-s).max());energy_change=None if last_energy is None else abs(met['total_J']-last_energy)/max(abs(last_energy),c['mechanical_equilibrium']['energy_floor_J'])
            phase_energy=met['total_J']-answer['before']['total_J']
            history.append(dict(step=step,strain=float(strain),iteration=iteration,phase_increment=ds,
                relative_energy_change=energy_change,phase_energy_change_J=phase_energy,**met))
            if phase_energy>c['tolerances']['phase_energy_absolute_J']+c['tolerances']['phase_energy_relative']*abs(answer['before']['total_J']):
                status='STOPPED';reason='PHASE_ENERGY_AUDIT_FAILED';break
            s=new.copy();last_energy=met['total_J']
            if (ds<=c['tolerances']['phase_step'] and energy_change is not None and
                energy_change<=c['tolerances']['stagger_energy_relative'] and
                met['force_residual']<=c['mechanical_equilibrium']['relative_force_residual_tolerance']):
                converged=True;accepted.append(dict(strain=float(strain),iterations=iteration+1,**met));break
        if status=='STOPPED':break
        if not converged:status='STOPPED';reason='STAGGERED_MAXIMUM_ITERATIONS_EXCEEDED';break
    np.savez_compressed(folder/'FINAL_STATE.npz',x=x,s=s,lower=prev['s_pf_bond'])
    result=dict(substeps=count,label='STEP_SUBDIVISION_DIAGNOSTIC_ONLY',status=status,reason=reason,
                endpoint_reached=bool(len(accepted)==count),accepted_substeps=accepted,
                stopped_or_last_strain=float(strain),last_computed_metrics=history[-1] if history else None)
    save_json(folder/'SUMMARY.json',result);save_json(folder/'HISTORY.json',history)
    return result


def main():
    a=audit_config();budget=worker_budget();save_json(OUT/'provenance/ENERGY_WORKER_BUDGET.json',budget)
    with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
        gradients=list(pool.map(gradient_audit,['A','B','failed_load_initial']))
    passed=all(g['passed'] for g in gradients)
    save_json(OUT/'monolithic/GRADIENT_GATE.json',dict(passed=passed,cases=gradients))
    print('GRADIENT_GATE',passed,gradients,flush=True)
    if not passed:raise SystemExit('MONOLITHIC_GRADIENT_NOT_VERIFIED')
    tasks=[]
    for key in ['A','B','C']:
        p,data=load_audit_state(key)
        tasks.extend((key,data['x'],data['s'],p.strain,p.lower,e) for e in a['hessian_epsilon'])
    with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
        Hessians=list(pool.map(hessian_audit,tasks))
    save_json(OUT/'hessian/BASE_STATES.json',Hessians)
    print('BASE HESSIANS',[(h['state'],h['epsilon'],h['minimum_reduced_eigenvalue']) for h in Hessians],flush=True)
    p,starts=deterministic_starts()
    jobs=[(method,label,x,s,p.strain,p.lower) for method in a['optimizer_methods'] for label,(x,s) in starts.items()]
    for key in ['A','B']:
        problem,data=load_audit_state(key)
        jobs.append(('L-BFGS-B',f'control_{key}',data['x'],data['s'],problem.strain,problem.lower))
    with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
        optimizers=list(pool.map(optimize_case,jobs))
    save_json(OUT/'monolithic/ALL_OPTIMIZERS.json',optimizers)
    print('OPTIMIZERS',[(r['method'],r['start'],r['optimizer_success'],r['message']) for r in optimizers],flush=True)
    ht=[]
    for row in optimizers:
        if row['first_order_candidate_passed']:
            state=read_npz(OUT/'monolithic'/row['method']/row['start']/'STATE.npz')
            ht.extend((row['method']+'_'+row['start'],state['x'],state['s'],row['strain'],state['lower'],e) for e in a['hessian_epsilon'])
    with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
        optimized_hessians=list(pool.map(hessian_audit,ht))
    save_json(OUT/'hessian/OPTIMIZER_STATES.json',optimized_hessians)
    with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
        substeps=list(pool.map(subdivision_case,a['subdivisions']))
    save_json(OUT/'continuation/ALL_SUBSTEPS.json',substeps)
    print('SUBSTEPS',[(r['substeps'],r['status'],r['reason']) for r in substeps],flush=True)


if __name__=='__main__':main()
