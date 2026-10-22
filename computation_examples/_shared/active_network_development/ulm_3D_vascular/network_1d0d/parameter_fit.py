"""Pure-0D weighted least squares with explicit observational rank guards.

Only geometry-network states enter residuals. Neither CFD fields, execution
modules, surrogates nor particle trajectories are accepted target sources.
"""
from __future__ import annotations
from dataclasses import dataclass, replace, asdict
from pathlib import Path
import csv
import json
import numpy as np
from scipy.optimize import least_squares

from .hydraulic_resistance import ROI_TARGET_Q_M3_S
from .parameterized_hydraulics import (
    MODEL_NAME, PORT_ORDER, OUTLET_ORDER, ParameterizedHydraulicSpec,
    solve_parameterized_operating_point)

PARAMETER_ORDER = ('s_O1','s_O2','terminal_resistance_O3_pa_s_m3')
OBSERVATION_KINDS = ('flow_fraction','outlet_flow_m3_s','real_cut_pressure_pa','pressure_difference_pa')


class UnidentifiableParameterizationError(ValueError):
    def __init__(self,audit):
        self.audit=audit
        super().__init__('Unidentifiable parameterization: '+json.dumps(audit,sort_keys=True))


@dataclass(frozen=True)
class Observation:
    kind: str
    port: str
    value: float
    sigma: float
    other_port: str | None = None

    def validate(self):
        if self.kind not in OBSERVATION_KINDS or self.port not in PORT_ORDER:
            raise ValueError('Unsupported observation kind or port')
        if not np.isfinite(self.value) or not np.isfinite(self.sigma) or self.sigma<=0:
            raise ValueError('Observation value finite and measurement sigma positive required')
        if self.kind in ('flow_fraction','outlet_flow_m3_s') and self.port not in OUTLET_ORDER:
            raise ValueError('Flow targets must be outlet observations; inlet Q is prescribed')
        if self.kind=='pressure_difference_pa':
            if self.other_port not in PORT_ORDER or self.port==self.other_port:
                raise ValueError('Pressure difference requires two distinct named real-cut ports')
        elif self.other_port is not None:
            raise ValueError('other_port only valid for pressure difference')

    @property
    def label(self):
        prefix={'flow_fraction':'flow_residual','outlet_flow_m3_s':'flow_rate_residual',
                'real_cut_pressure_pa':'pressure_residual','pressure_difference_pa':'pressure_difference_residual'}[self.kind]
        return prefix+'_'+self.port+('_minus_'+self.other_port if self.other_port else '')+'_dimensionless'

    def predicted(self,result):
        i=PORT_ORDER.index(self.port)
        if self.kind=='flow_fraction': return result.signed_outlet_fractions[i-1]
        if self.kind=='outlet_flow_m3_s': return result.port_flow_m3_s[i]
        if self.kind=='real_cut_pressure_pa': return result.port_pressure_pa[i]
        return result.port_pressure_pa[i]-result.port_pressure_pa[PORT_ORDER.index(self.other_port)]


@dataclass(frozen=True)
class Targets:
    kind: str
    observations: tuple[Observation,...]
    evidence_reference: str

    def validate(self):
        if self.kind not in ('SYNTHETIC_TARGET','EXPERIMENTAL_TARGET'):
            raise ValueError('Only SYNTHETIC_TARGET from pure 0D or EXPERIMENTAL_TARGET from measurements is allowed; CFD/particle-derived targets forbidden')
        if not self.evidence_reference or not self.observations:
            raise ValueError('Explicit target evidence and nonempty observations required')
        for o in self.observations:o.validate()
        if len(set(o.label for o in self.observations))!=len(self.observations):
            raise ValueError('Duplicate observation')

    def canonical(self):
        self.validate()
        return replace(self,observations=tuple(sorted(self.observations,key=lambda o:(
            OBSERVATION_KINDS.index(o.kind),PORT_ORDER.index(o.port),o.other_port or ''))))


@dataclass(frozen=True)
class LogPrior:
    parameter: str
    mean_log_dimensionless: float
    sigma_log_dimensionless: float


def _structural_identifiability(targets,free,priors,spec):
    # At fixed Qin, Q targets and fractions constrain the SAME two-dimensional
    # outlet-flow state. Counting all Q and fraction rows independently is wrong.
    flow_rows=[]; pressure_rows=[]
    basis={'O1':[1.,0.],'O2':[0.,1.],'O3':[-1.,-1.]}
    for obs in targets.observations:
        if obs.kind in ('flow_fraction','outlet_flow_m3_s'):flow_rows.append(basis[obs.port])
        else:
            row=np.zeros(4);row[PORT_ORDER.index(obs.port)]=1
            if obs.other_port:row[PORT_ORDER.index(obs.other_port)]-=1
            if spec.o3_mode=='legacy_reference':row[3]=0 # O3 pressure is fixed, not informative
            pressure_rows.append(row)
    rank=lambda a:int(np.linalg.matrix_rank(a)) if len(a) else 0
    flow_rank, pressure_rank=rank(flow_rows),rank(pressure_rows)
    prior_rank=len(priors)
    audit=dict(number_of_free_parameters=len(free),free_parameters=list(free),
               number_of_independent_flow_constraints=flow_rank,pressure_constraints=pressure_rank,
               prior_constraints=prior_rank,
               fixed_anchor_parameters=[p for p in PARAMETER_ORDER if p not in free
                    and (p!='terminal_resistance_O3_pa_s_m3' or spec.o3_mode=='terminal_resistance')],
               inactive_parameters=['terminal_resistance_O3_pa_s_m3'] if spec.o3_mode=='legacy_reference' else [],
               observed_flow_rows=len(flow_rows),observed_pressure_rows=len(pressure_rows),
               data_constraints_upper_bound=flow_rank+pressure_rank,
               qualification='count is only necessary; finite-difference Jacobian rank is also checked')
    if len(free)>flow_rank+pressure_rank+prior_rank:
        raise UnidentifiableParameterizationError(audit|{'reason':'free parameters exceed independent observations plus explicit log priors'})
    return audit


def _rank_audit(jac):
    singular=np.linalg.svd(jac,compute_uv=False)
    tol=max(1e-8,float(singular[0])*1e-7) if len(singular) else 1e-8
    return dict(rank=int(np.sum(singular>tol)),singular_values=singular.tolist(),
                normalized_residual_jacobian_rank_tolerance=tol)


def fit_parameters(cache,initial_spec,targets,*,free_parameters=('s_O1','s_O2'),
                   target_roi_flow_m3_s=ROI_TARGET_Q_M3_S,priors=(),
                   R_ref_pa_s_m3=1e17,physical_bounds=None,max_nfev=1000,
                   trace_path=None,provenance=None,state_output_path=None):
    initial_spec.validate();targets=targets.canonical()
    if len(set(free_parameters))!=len(free_parameters) or any(p not in PARAMETER_ORDER for p in free_parameters):
        raise ValueError('Unknown or duplicate free parameter')
    free=tuple(p for p in PARAMETER_ORDER if p in free_parameters)
    if not free: raise ValueError('Use the forward solver for a model with no free parameters')
    if 'terminal_resistance_O3_pa_s_m3' in free and initial_spec.o3_mode!='terminal_resistance':
        raise ValueError('O3 resistance is free only in terminal_resistance mode')
    if not np.isfinite(R_ref_pa_s_m3) or R_ref_pa_s_m3<=0:
        raise ValueError('R_ref_pa_s_m3 must be finite positive')
    if not isinstance(max_nfev,int) or isinstance(max_nfev,bool) or max_nfev<=0:
        raise ValueError('max_nfev must be a positive integer')
    priors=tuple(sorted(priors,key=lambda p:PARAMETER_ORDER.index(p.parameter) if p.parameter in PARAMETER_ORDER else -1))
    if len(set(p.parameter for p in priors))!=len(priors):raise ValueError('Duplicate prior')
    for p in priors:
        if (p.parameter not in free or not np.isfinite(p.mean_log_dimensionless)
                or not np.isfinite(p.sigma_log_dimensionless) or p.sigma_log_dimensionless<=0):
            raise ValueError('Prior must constrain a free parameter with finite mean and positive sigma')
    audit=_structural_identifiability(targets,free,priors,initial_spec)
    references=np.array([R_ref_pa_s_m3 if p=='terminal_resistance_O3_pa_s_m3' else 1. for p in free])
    x0=np.log(np.array([getattr(initial_spec,p) for p in free])/references)
    lower=np.full(len(free),-np.inf);upper=-lower
    for parameter,bounds in (physical_bounds or {}).items():
        if parameter not in free:raise ValueError('Bounds must refer to a free parameter')
        lo,hi=bounds
        if not np.isfinite([lo,hi]).all() or not 0<lo<hi:raise ValueError('Bounds must be finite, positive and ordered')
        i=free.index(parameter);lower[i],upper[i]=np.log(np.array([lo,hi])/references[i])
    if np.any(x0<lower) or np.any(x0>upper):raise ValueError('Initial parameter lies outside physical bounds')
    trace=[];phase='initial'
    residual_names=[o.label for o in targets.observations]+['prior_residual_'+p.parameter+'_dimensionless' for p in priors]

    def decode(x):
        with np.errstate(over='raise',invalid='raise',under='raise'):
            try:values=references*np.exp(x)
            except FloatingPointError as exc:raise ValueError('Log parameters exceeded finite positive physical range') from exc
        return replace(initial_spec,**dict(zip(free,values)))

    def evaluate(x):
        row=dict(evaluation=len(trace)+1,phase=phase,status='FAILED',error='')
        try:
            spec=decode(x)
            row.update(s_O1_dimensionless=spec.s_O1,s_O2_dimensionless=spec.s_O2,
                       R_terminal_O3_pa_s_m3=spec.terminal_resistance_O3_pa_s_m3,mu_pa_s=spec.mu_pa_s)
            result=solve_parameterized_operating_point(cache,spec,target_roi_flow_m3_s)
            data=[(o.predicted(result)-o.value)/o.sigma for o in targets.observations]
            prior=[(x[free.index(p.parameter)]-p.mean_log_dimensionless)/p.sigma_log_dimensionless for p in priors]
            residual=np.array(data+prior)
            if not np.all(np.isfinite(residual)):raise ValueError('Nonfinite residual')
            row.update(status='PASS',mass_residual_dimensionless=result.mass_audit['max_relative_residual'],
                       objective_norm_dimensionless=float(np.linalg.norm(residual)))
            for p,q,pressure in zip(PORT_ORDER,result.port_flow_m3_s,result.port_pressure_pa):
                row['Q_'+p+'_m3_s']=float(q);row['P_'+p+'_pa']=float(pressure)
            for p,f in zip(OUTLET_ORDER,result.signed_outlet_fractions):row['f_'+p+'_dimensionless']=float(f)
            row.update(dict(zip(residual_names,residual.tolist())))
            trace.append(row)
            return residual,result
        except Exception as exc:
            row['error']=type(exc).__name__+': '+str(exc);trace.append(row)
            raise

    def jacobian(x,label):
        nonlocal phase
        phase=label;columns=[]
        # Central finite difference of actual whole-network evaluations. No CFD
        # and no analytical outlet black-box replaces the forward model.
        for j in range(len(free)):
            step=1e-4*max(1.,abs(x[j]));lo=x.copy();hi=x.copy()
            lo[j]=max(lower[j],x[j]-step);hi[j]=min(upper[j],x[j]+step)
            columns.append((evaluate(hi)[0]-evaluate(lo)[0])/(hi[j]-lo[j]))
        return np.column_stack(columns)

    try:
        initial_residual,_=evaluate(x0)
        j0=jacobian(x0,'initial_identifiability')
        audit['initial_data_jacobian']=_rank_audit(j0[:len(targets.observations)])
        audit['initial_augmented_jacobian']=_rank_audit(j0)
        if audit['initial_augmented_jacobian']['rank']<len(free):
            raise UnidentifiableParameterizationError(audit|{'reason':'initial normalized Jacobian is numerically rank deficient'})
        phase='optimizer'
        opt=least_squares(lambda x:evaluate(x)[0],x0,jac='3-point',bounds=(lower,upper),
                          max_nfev=max_nfev,ftol=1e-11,xtol=1e-11,gtol=1e-11)
        jf=jacobian(opt.x,'final_identifiability')
        audit['final_data_jacobian']=_rank_audit(jf[:len(targets.observations)])
        audit['final_augmented_jacobian']=_rank_audit(jf)
        audit['interpretation']=('LOCALLY_FULL_RANK_DATA_JACOBIAN_NOT_GLOBAL_UNIQUENESS'
            if audit['final_data_jacobian']['rank']==len(free) else 'PRIOR_REGULARIZED_NOT_DATA_IDENTIFIED')
        phase='final';final_residual,final=evaluate(opt.x)
        if audit['final_augmented_jacobian']['rank']<len(free):
            raise UnidentifiableParameterizationError(audit|{'reason':'final augmented Jacobian is rank deficient'})
        if state_output_path is not None:
            np.savez_compressed(state_output_path,original_node_ids=cache.node_ids,pressure_pa=final.pressure_pa,
                edge_nodes=final.edge_nodes,edge_flow_m3_s=final.edge_flow_m3_s,node_outflow_m3_s=final.node_outflow_m3_s,
                edge_resistance_pa_s_m3=final.edge_resistance_pa_s_m3,original_node_radius_m=final.effective_radius_m,
                port_indices=cache.port_indices,port_names=np.array(PORT_ORDER),
                augmented_network_metadata_json=json.dumps(final.closure_provenance))
        return dict(status='FIT_CONVERGED' if opt.success else 'FIT_NOT_CONVERGED',model_name=MODEL_NAME,
            parameterization={p:'log(R/R_ref)' if p.startswith('terminal') else 'log(radius_scale)' for p in free},
            R_ref_pa_s_m3=R_ref_pa_s_m3,parameter_source=targets.kind,
            initial_parameters=asdict(initial_spec),final_parameters=asdict(final.spec),
            roi_target_flow_m3_s=target_roi_flow_m3_s,
            targets=asdict(targets),priors=[asdict(p) for p in priors],
            initial_residual=initial_residual.tolist(),final_residual=final_residual.tolist(),
            residual_labels=residual_names,initial_objective_norm=float(np.linalg.norm(initial_residual)),
            final_objective_norm=float(np.linalg.norm(final_residual)),
            optimizer_status=dict(success=bool(opt.success),status=int(opt.status),message=opt.message,
                                  optimality=float(opt.optimality),active_mask=opt.active_mask.tolist()),
            nfev=int(opt.nfev),njev=int(opt.njev),total_forward_evaluations=len(trace),
            identifiability_audit=audit,input_hashes=dict(cache.input_hashes),
            source_geometry_sha256=cache.geometry_sha256,provenance=dict(provenance or {}),
            prediction=final.summary(),CFD_calls=0,particle_simulation_calls=0,
            scientific_limit='targets must be genuinely synthetic 0D or independent measurements; provenance declaration cannot authenticate external measurements')
    finally:
        if trace_path is not None:
            fields=list(dict.fromkeys(k for row in trace for k in row))
            with Path(trace_path).open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(trace)
