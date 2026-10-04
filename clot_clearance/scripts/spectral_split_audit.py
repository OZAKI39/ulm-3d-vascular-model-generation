"""Spectral candidate audit; all existing numerical/model files remain untouched."""
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
import ast,copy,json,hashlib,time,os,subprocess,resource
import numpy as np
from scipy.linalg import sqrtm
from scripts import spatial_audit_v2 as v2
from pd_clot.pfpd_candidate import baqp
from pd_clot.pfpd_candidate import branch_stability_tracking as tracker
from pd_clot.pfpd_candidate import coupling_stability_audit as ca
from pd_clot.pfpd_candidate.svk_spectral_split import (green_lagrange,spectral_green_split,
    svk_split_energy,svk_split_second_pk,svk_split_response,degraded_svk_split_response)
from pd_clot.pfpd_candidate.energy import Material,constitutive
from pd_clot.pfpd_candidate.at1_evolution import at1_constrained_update

ROOT=v2.ROOT;OUT=ROOT/'verification/pfpd_at1_spectral_split_v1';RESULTS=ROOT/'results/pfpd_at1_spectral_split_v1'
read=ca.read_json;save=ca.save_json
STRAINS=[.00025,.005,.009,.01]
WINDOW=[.009,.0095,.00975,.009875,.0099375,.00996875,.009984375,.0099921875,.01,.010125,.01025,.0105]

def declare():
    cfg=read(v2.OUT/'provenance/PREDECLARED_CONFIG.json')
    cfg.update(task='PFPD_AT1_SPECTRAL_SPLIT_V1',intact_strains=STRAINS,branch_window=WINDOW,
        decomposition_relative_tolerance=2e-12,stress_decomposition_relative_tolerance=2e-12,
        energy_derivative_relative_tolerance=5e-6,energy_derivative_epsilon=[1e-5,3e-6,1e-6],
        intact_equivalence_relative_tolerance=2e-12,tensor_oracle_relative_tolerance=2e-8,
        finite_J_minimum=0.,pure_tension_policy='DIAGNOSTIC_ONLY_NEW_CANDIDATE; old tolerance and old gate unchanged',
        negative_fraction_definition='Fraction with minimum principal Green strain < -1e-9 (historical threshold for diagnostic counting only); raw <0 also saved.',
        physical_constants_source=str(ROOT/'verification/pfpd_at1_evolution_v1/provenance/PREDECLARED_CONFIG.json'),
        cpu_backend='NumPy batched eigh; existing NumPy BAQP; existing ADR and branch tracker unchanged',
        gpu_backend='NOT_USED_NO_EXISTING_BATCH_LINALG_BACKEND',gpu_benchmark='NOT_RUN',
        GPU_note='RTX 4060 and numba.cuda available, but CuPy/JAX/Torch unavailable. No custom GPU eigensolver or CUDA stack installation for this audit.',
        intact_solver='Existing adr_adapter_at1.equilibrate_tensile; fixed s=0 and original numerical settings. Intact identity established before old SVK force callback is reused.',
        gates=['spectral constitutive','s=0 equivalence','intact equilibrium C3/C5/C7','spectral AT1 branch','spatial stability'])
    save(OUT/'provenance/PREDECLARED_CONFIG.json',cfg)

def tensor_cases():
    return dict(pure_tension=np.diag([.01,.02,.03]),pure_compression=np.diag([-.01,-.02,-.03]),
        mixed=np.diag([-.02,.01,.03]),shear=np.array([[0,.02,0],[.02,0,0],[0,0,0]]),
        hydrostatic_tension=np.eye(3)*.02,hydrostatic_compression=-np.eye(3)*.02,
        repeated_eigenvalues=np.diag([.02,.02,-.01]),near_zero=np.diag([-1e-13,.01,.02]))

def F_from_E(E):
    e,Q=np.linalg.eigh(np.eye(3)+2*E)
    if np.any(e<=0):raise ValueError('Non-admissible prescribed strain')
    return (Q*np.sqrt(e)[...,None,:])@Q.swapaxes(-1,-2)

def oracle_functions():
    path=ROOT/'verification/pfpd_at1_evolution_v1/resources/external/phasefieldx/src/phasefieldx/Element/Phase_Field_Fracture/split_energy_stress_tangent_functions.py'
    names=['psi_a_spectral','psi_b_spectral','sigma_a_spectral','psi_a_deviatoric','psi_b_deviatoric']
    tree=ast.parse(path.read_text());nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    def absolute(E):return np.real_if_close(sqrtm(E@E)).real
    namespace=dict(ufl=SimpleNamespace(tr=np.trace,inner=lambda a,b:np.sum(a*b),Identity=np.eye),epsilon=lambda E:E,
        macaulay_bracket_positive=lambda a:.5*(a+abs(a)),macaulay_bracket_negative=lambda a:.5*(a-abs(a)),
        spectral_positive_part=lambda E:.5*(E+absolute(E)),spectral_negative_part=lambda E:.5*(E-absolute(E)),
        deviatoric_part=lambda E:E-np.eye(3)*np.trace(E)/3)
    # Execute only the inspected small algebraic functions, not package imports or FEM solver.
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),namespace)
    return namespace,path

def verify_constitutive():
    cfg=read(OUT/'provenance/PREDECLARED_CONFIG.json');rng=np.random.default_rng(6093);rows=[];derivatives=[];oracle=[]
    op,source=oracle_functions()
    for mat in [Material('svk',1000.,666.6666666666666),Material('svk',1000.,2300.)]:
        for name,E in tensor_cases().items():
            F=F_from_E(E);r=svk_split_response(F,mat);W,P=constitutive(F,mat);S=mat.lame*np.trace(r.E)*np.eye(3)+2*mat.mu*r.E
            row=dict(case=name,lambda_Pa=mat.lame,Psi_plus=float(r.Psi_plus),Psi_minus=float(r.Psi_minus),
                energy_decomposition_error=float(abs(r.Psi_total-W)/max(abs(W),1e-14)),
                stress_decomposition_error=float(np.linalg.norm(r.S_total-S)/max(np.linalg.norm(S),1e-14)),
                F_PK_error=float(np.linalg.norm(r.P_total-P)/max(np.linalg.norm(P),1e-14)))
            rows.append(row)
            pp=op['psi_a_spectral'](E,mat.lame,mat.mu);pm=op['psi_b_spectral'](E,mat.lame,mat.mu);sp=op['sigma_a_spectral'](E,mat.lame,mat.mu)
            cp,cm=svk_split_energy(E,mat.lame,mat.mu);cs,_=svk_split_second_pk(E,mat.lame,mat.mu)
            oracle.append(dict(case=name,lambda_Pa=mat.lame,Psi_plus_error=float(abs(cp-pp)/max(abs(cp),1)),
                Psi_minus_error=float(abs(cm-pm)/max(abs(cm),1)),S_plus_error=float(np.linalg.norm(cs-sp)/max(np.linalg.norm(cs),1)),
                spectral_plus=float(cp),spectral_minus=float(cm),Amor_plus=float(op['psi_a_deviatoric'](E,mat.lame,mat.mu)),Amor_minus=float(op['psi_b_deviatoric'](E,mat.lame,mat.mu))))
            if np.min(abs(np.linalg.eigvalsh(E)))>1e-4 and abs(np.trace(E))>1e-4:
                for k in range(3):
                    H=rng.normal(size=(3,3));H=(H+H.T)/2;H/=np.linalg.norm(H)
                    for s in [0.,.37]:
                        Sp,Sm=svk_split_second_pk(E,mat.lame,mat.mu);analytic=float(np.sum(((1-s)**2*Sp+Sm)*H))
                        for eps in cfg['energy_derivative_epsilon']:
                            def energy(T):
                                a,b=svk_split_energy(T,mat.lame,mat.mu);return (1-s)**2*a+b
                            fd=float((energy(E+eps*H)-energy(E-eps*H))/(2*eps))
                            derivatives.append(dict(case=name,lambda_Pa=mat.lame,s=s,direction=k,epsilon=eps,analytic=analytic,fd=fd,
                                relative_error=abs(fd-analytic)/max(abs(analytic),np.linalg.norm((1-s)**2*Sp+Sm),1e-12)))
    # One-sided switching diagnostic: energy is C1; its tangent need not be smooth.
    E=np.diag([0.,.01,-.02]);H=np.diag([1.,0,0]);switch=[];mat=Material('svk',1000.,2300.)
    for sign in [-1,1]:
        for eps in cfg['energy_derivative_epsilon']:
            f=lambda T:sum(svk_split_energy(T,mat.lame,mat.mu))
            switch.append(dict(sign=sign,epsilon=eps,directional_derivative=float((f(E+sign*eps*H)-f(E))/(sign*eps))))
    invariance=[]
    for k in range(128):
        F=np.eye(3)+.13*rng.normal(size=(3,3))
        if np.linalg.det(F)<=0:raise RuntimeError('Deterministic F fixture unexpectedly inadmissible')
        mat=Material('svk',1000.,666.6666666666666 if k%2 else 2300.);r=degraded_svk_split_response(F,0.,mat);W,P=constitutive(F,mat)
        S=np.linalg.solve(F,P)
        invariance.append(dict(index=k,energy=float(abs(r.Psi-W)/max(abs(W),1e-14)),
            stress=float(np.linalg.norm(r.S-S)/max(np.linalg.norm(S),1e-14)),PK=float(np.linalg.norm(r.P-P)/max(np.linalg.norm(P),1e-14))))
    maxima=dict(Psi_decomposition=max(r['energy_decomposition_error'] for r in rows),stress_decomposition=max(r['stress_decomposition_error'] for r in rows),
        energy_derivative=max(r['relative_error'] for r in derivatives),intact_equivalence=max(max(r[k] for k in ['energy','stress','PK']) for r in invariance),
        oracle=max(max(r[k] for k in ['Psi_plus_error','Psi_minus_error','S_plus_error']) for r in oracle))
    gate1=maxima['Psi_decomposition']<2e-12 and maxima['stress_decomposition']<2e-12 and maxima['energy_derivative']<5e-6 and maxima['oracle']<2e-8
    gate2=maxima['intact_equivalence']<2e-12
    result=dict(gate1_passed=gate1,gate2_passed=gate2,maxima=maxima,tensors=rows,derivatives=derivatives,switching_diagnostic=switch,intact_random_states=invariance)
    save(OUT/'constitutive/VERIFICATION.json',result)
    save(OUT/'references/PHASEFIELDX_TENSOR_ORACLE.json',dict(status='ORACLE_USED',commit='a9714c122497f8829860efc68d776e4324d48eca',source=str(source),
        sha256=hashlib.sha256(source.read_bytes()).hexdigest(),functions_executed=['psi_a_spectral','psi_b_spectral','sigma_a_spectral','psi_a_deviatoric','psi_b_deviatoric'],
        spectral_input_backend='Independent scipy.linalg.sqrtm(E@E), E+/- = (E +/- |E|)/2; numerical shim for UFL scalar/tensor operations; epsilon(u) replaced by supplied E',
        scope='Algebraic tensor-level comparison only; NOT finite-strain solver equivalence; no PhaseFieldX FEM solver imported or vendored',cases=oracle))
    print('CONSTITUTIVE',maxima,'GATE1',gate1,'GATE2',gate2,flush=True)
    assert gate1 and gate2
    return result


def split_baqp(family,x,material,s):
    """Reuse the intact BAQP assembly with a scoped constitutive callback.

    Mechanical energy/PK are already degraded in the callback, so outer g=1.
    W0 is the compatibility alias for active driver Psi_plus, never Psi_total.
    The scope is serial/process-local and restores the old callback on exit.
    """
    responses=[]
    def callback(F,mat):
        r=degraded_svk_split_response(F,s,mat);responses.append(r);return r.Psi,r.P
    with patch.object(baqp,'constitutive',callback):
        state=baqp.evaluate(family,x,material,g=np.ones(len(family.pairs)),h=np.ones(len(family.pairs)))
    r=responses[0];state.W0=r.Y;state.split_response=r
    return state


class SpectralProblem(v2.Problem):
    def state(self,x,s):
        if hasattr(self,'_last_split') and np.array_equal(x,self._last_x) and np.array_equal(s,self._last_s):return self._last_split
        state=split_baqp(self.family,x,self.material,s)
        self._last_x=np.array(x,copy=True);self._last_s=np.array(s,copy=True);self._last_split=state
        return state

    def metrics(self,x,s):
        result=super().metrics(x,s);state=self.state(x,s);sp=state.split_response.split
        positive=float(self.weights@sp.Psi_plus);negative=float(self.weights@sp.Psi_minus);den=positive+negative
        result.update(min_principal_green_strain=float(sp.eigenvalues.min()),
            negative_principal_fraction=float(np.mean(sp.eigenvalues.min(axis=1)<-1e-9)),
            total_Psi_plus_J=positive,total_Psi_minus_J=negative,
            Psi_plus_fraction=positive/den if den else None,Psi_minus_fraction=negative/den if den else None)
        return result

def spectral_stationary(p,x,s):
    """Scoped domain adapter around the UNCHANGED verified branch solver.

    The old residual/KKT/Jacobian algorithms and all tolerances are retained.
    Only the new constitutive candidate's Green-strain rejection is removed.
    """
    cfg=read(OUT/'provenance/PREDECLARED_CONFIG.json');original=tracker.stationary_residual;last_log=[0.]
    if p.n>3:cfg['least_squares']['tr_solver']='lsmr'
    def mixed_residual(problem,z,check_domain=True):
        if check_domain:
            try:
                xx,ss=problem.unpack(z);m=problem.metrics(xx,ss)
            except (ValueError,FloatingPointError) as exc:
                raise tracker.BranchDomainExit('FINITE_OR_INVERSION_SAFETY_GATE',z,dict(error=str(exc))) from exc
            if m['maximum_s']>=problem.sc:raise tracker.BranchDomainExit('HIGH_DAMAGE_GATE_REACHED',z,m)
            if m['minimum_lower_margin']<-cfg['bounds_tolerance']:raise tracker.BranchDomainExit('IRREVERSIBILITY_BOUND_VIOLATION',z,m)
            if time.time()-last_log[0]>30:
                last_log[0]=time.time();print(f'BRANCH C{p.n} eps={p.strain:g} force={m["force_residual"]:.4g} KKT={m["projected_phase_KKT"]:.4g} minE={m["minimum_Green_strain"]:.5g}',flush=True)
        return original(problem,z,check_domain=False)
    with patch.object(tracker,'settings',lambda:cfg),patch.object(tracker,'stationary_residual',mixed_residual):
        if p.n>3:
            with patch.object(tracker,'residual_jacobian',v2.matrix_free_residual_jacobian):answer=tracker.solve_stationary_branch_state(p,x,s)
        else:answer=tracker.solve_stationary_branch_state(p,x,s)
    answer['old_pure_tension_gate_diagnostic']=answer['gates'].pop('pure_tension')
    state=p.state(answer['x'],answer['s']);J=np.linalg.det(state.Ftilde)
    answer['gates']['positive_finite_J']=bool(np.isfinite(J).all() and np.all(J>0))
    answer['accepted']=bool(all(answer['gates'].values()));answer['candidate']='PFPD_AT1_SPECTRAL_SPLIT_V1'
    return answer

def mixed_probe_point(p,z,direction,amplitude):
    """Same V1 cone probe with only the authorized physical-domain change."""
    trial=z+amplitude*direction;x,s=p.unpack(trial)
    if not (np.all(s>=p.lower) and np.all(s<p.sc)):
        return dict(amplitude=amplitude,feasible=False,reason='PHASE_BOUNDS')
    try:state=p.state(x,s)
    except (ValueError,FloatingPointError) as exc:
        return dict(amplitude=amplitude,feasible=False,reason='FINITE_OR_INVERSION_SAFETY_GATE',detail=str(exc))
    J=np.linalg.det(state.Ftilde);finite=bool(np.isfinite(J).all() and np.all(J>0))
    fixed=bool(np.array_equal(x[p.cloud.fixed],p.boundary[p.cloud.fixed]))
    return dict(amplitude=amplitude,feasible=finite and fixed,minimum_Green_strain=ca.green_minimum(state),
        minimum_J=float(J.min()),fixed_boundary=fixed,reason=None if finite and fixed else 'FINITE_OR_FIXED_BOUNDARY',
        green_strain_policy='DIAGNOSTIC_ONLY')

def constitutive_baqp_fixtures():
    p=SpectralProblem(3,0);rows=[]
    for name,F in [('compression',np.eye(3)*.99),('shear',np.array([[1,.02,0],[0,1,0],[0,0,1]]))]:
        x=p.family.X@F.T;state=p.state(x,p.lower);r=state.split_response.split
        update=at1_constrained_update(p.lower,state.W0,p.yc,p.sc)
        row=dict(case=name,Psi_plus_min=float(r.Psi_plus.min()),Psi_minus_min=float(r.Psi_minus.min()),
            min_principal_green_strain=float(r.eigenvalues.min()),driver_max=float(state.W0.max()),
            maximum_phase_increase=float((update.s_new-p.lower).max()),mixed_state_admitted=True)
        rows.append(row)
    assert rows[0]['driver_max']<1e-20 and rows[0]['maximum_phase_increase']==0
    assert rows[1]['Psi_plus_min']>0 and rows[1]['Psi_minus_min']>0
    # Reuse saved pure-tension endpoints only for cheap constitutive-limit checks.
    regression=[]
    for e in [.0099921875,.01]:
        p=SpectralProblem(3,e);x=v2.affine(p);old=baqp.evaluate(p.family,x,p.material);new=p.state(x,p.lower)
        regression.append(dict(strain=e,energy_relative_error=abs(new.energy-old.energy)/old.energy,
            reaction_relative_error=float(abs(new.force[p.top,2].sum()-old.force[p.top,2].sum())/abs(old.force[p.top,2].sum())),
            driver_relative_error=float(np.max(abs(new.W0-old.W0))/max(old.W0.max(),1e-30)),
            Psi_minus_max=float(new.split_response.split.Psi_minus.max())))
    save(OUT/'constitutive/BAQP_SMALL_FIXTURES.json',dict(fixtures=rows,C3_pure_tension_limit=regression))
    assert max(max(r[k] for k in ['energy_relative_error','reaction_relative_error','driver_relative_error']) for r in regression)<2e-12

def spatial_diagnostics(p,x,state):
    split=svk_split_response(state.Ftilde,p.material);point=svk_split_response(state.Fbar,p.material)
    eb=split.eigenvalues.min(axis=1);ep=point.eigenvalues.min(axis=1);weights=p.weights
    pos=float(weights@split.Psi_plus);neg=float(weights@split.Psi_minus);den=pos+neg
    mid=p.family.X[p.family.pairs].mean(axis=1);bi=int(np.argmin(eb));pi=int(np.argmin(ep));X=p.family.X
    def location(a):return dict(position_m=a.tolist(),clamp_distance_m=float(v2.L/2-abs(a[2])),side_distance_m=float(v2.L/2-max(abs(a[:2]))))
    groups={};h=p.cloud.spacing
    for name,mask in dict(near_top=X[:,2]>=v2.L/2-h/2,near_bottom=X[:,2]<=-v2.L/2+h/2,
            lateral=np.max(abs(X[:,:2]),axis=1)>=v2.L/2-h/2,interior_like=np.max(abs(X),axis=1)<v2.L/2-h/2).items():
        groups[name]=dict(count=int(mask.sum()),min=float(ep[mask].min()) if mask.any() else None,
            percentiles=np.percentile(ep[mask],[0,1,5,25,50,75,95,99,100]).tolist() if mask.any() else [],
            negative_fraction=float(np.mean(ep[mask]<-1e-9)) if mask.any() else None)
    row=dict(min_principal_green_strain=float(eb.min()),minimum_particle_green_strain=float(ep.min()),
        negative_principal_fraction=float(np.mean(eb<-1e-9)),negative_principal_particle_fraction=float(np.mean(ep<-1e-9)),
        negative_principal_bond_fraction=float(np.mean(eb<-1e-9)),negative_raw_particle_fraction=float(np.mean(ep<0)),negative_raw_bond_fraction=float(np.mean(eb<0)),
        total_Psi_plus_J=pos,total_Psi_minus_J=neg,Psi_plus_fraction=pos/den if den else None,Psi_minus_fraction=neg/den if den else None,
        most_compressive_particle=dict(ID=pi,**location(X[pi])),most_compressive_BAQP=dict(ID=bi,**location(mid[bi])),
        particle_percentiles=np.percentile(ep,[0,1,5,25,50,75,95,99,100]).tolist(),bond_percentiles=np.percentile(eb,[0,1,5,25,50,75,95,99,100]).tolist(),
        spatial_groups=groups,group_note='Groups overlap at clamp/lateral edges; fractions are within each group.')
    return row,dict(particle_min_green=ep,bond_min_green=eb,Psi_plus=split.Psi_plus,Psi_minus=split.Psi_minus,
        particle_X=X,bond_midpoints=mid)

def intact_equilibria():
    from pd_clot.pfpd_candidate import adr_adapter_at1 as adapter
    from pd_clot.regularization import equilibrium as mechanics
    from dataclasses import replace
    verification=read(OUT/'constitutive/VERIFICATION.json');assert verification['gate1_passed'] and verification['gate2_passed']
    cfg=ca.physics_config();allrows=[]
    for n in [3,5,7]:
        for e in STRAINS:
            folder=OUT/f'equilibrium/C{n}/{e:.8f}';folder.mkdir(parents=True,exist_ok=True)
            p=SpectralProblem(n,e);x0=v2.affine(p);t0=time.time();counter=[0];last=[0.]
            view=adapter.boundary_view(p.cloud,p.top,e)
            guard_cfg=dict(safety=dict(maximum_displacement_m=.01,maximum_radius_ratio=1.),pipe=dict(radius_m=1.))
            base=adapter.FrozenPhaseForce
            class LoggedForce(base):
                def __call__(self,x):
                    mechanics.check_positions(view,x,guard_cfg)
                    state=super().__call__(x);counter[0]+=1
                    if time.time()-last[0]>30:
                        last[0]=time.time();print(f'INTACT C{n} eps={e:g} force_calls={counter[0]} elapsed={time.time()-t0:.1f}s',flush=True)
                    return state
            print(f'START INTACT C{n} eps={e:g}',flush=True)
            with patch.object(adapter,'FrozenPhaseForce',LoggedForce):
                answer=adapter.equilibrate_tensile(p.cloud,p.family,p.material,p.top,x0,p.lower,e,cfg)
            x=answer.x;state=baqp.evaluate(p.family,x,p.material);diag,arrays=spatial_diagnostics(p,x,state)
            force=float(np.linalg.norm(state.force[p.free])/cfg['mechanical_equilibrium']['force_floor_N']);J=np.linalg.det(state.Ftilde)
            mechanics.check_positions(view,x,guard_cfg)
            gates=dict(existing_ADR_converged=answer.converged,force=force<=cfg['mechanical_equilibrium']['relative_force_residual_tolerance'],
                fixed_boundary=np.array_equal(x[p.cloud.fixed],view.X[p.cloud.fixed]),finite=bool(np.isfinite(x).all() and np.isfinite(J).all()),
                positive_J=bool(np.all(J>0)),safety=True,phase_frozen=bool(np.all(p.lower==0)))
            reaction=float(-state.force[p.top,2].sum())
            row=dict(level=f'C{n}',strain=e,accepted=bool(all(gates.values())),gates=gates,force_residual=force,reaction_N=reaction,nominal_stress_Pa=reaction/v2.L**2,
                elastic_energy_J=state.energy,minimum_J=float(J.min()),maximum_displacement_m=float(np.linalg.norm(x-p.family.X,axis=1).max()),
                maximum_deviation_from_affine_m=float(np.linalg.norm(x-x0,axis=1).max()),**diag,wall_seconds=time.time()-t0,
                force_evaluations=counter[0],peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,ADR=answer.summary)
            save(folder/'SUMMARY.json',row);save(folder/'ADR_HISTORY.json',answer.history)
            np.savez_compressed(folder/'STATE.npz',x=x,s=np.zeros(p.ns),x_affine=x0,force=state.force,**arrays)
            allrows.append(row);save(OUT/'equilibrium/ALL.json',allrows)
            print(f'END INTACT C{n} eps={e:g} accepted={row["accepted"]} force={force:.4g} minE={diag["min_principal_green_strain"]:.6g} wall={row["wall_seconds"]:.1f}s',flush=True)
            if not row['accepted']:
                save(OUT/'GATES.json',dict(GATE_1='PASSED',GATE_2='PASSED',GATE_3='FAILED',GATE_4='NOT_RUN',GATE_5='NOT_RUN',stop=f'C{n} strain {e}: intact equilibrium failed'))
                return allrows
    save(OUT/'GATES.json',dict(GATE_1='PASSED',GATE_2='PASSED',GATE_3='PASSED',GATE_4='PENDING',GATE_5='PENDING'))
    return allrows

if __name__=='__main__':
    declare();verify_constitutive();constitutive_baqp_fixtures()
