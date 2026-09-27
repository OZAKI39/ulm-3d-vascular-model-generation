"""Bounded pure-0D DESIGN optimization, separate from parameter identification.

P/Q are always full-network states. An active design bound is allowed; data
Jacobian rank is not an acceptance gate. No CFD/particle data or execution.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import os
import numpy as np
from scipy.optimize import least_squares
import scipy
import yaml

from .audit import sha256, write_json, write_csv
from .parameterized_hydraulics import (
    MODEL_NAME, PORT_ORDER, OUTLET_ORDER, ParameterizedHydraulicSpec,
    load_geometry_cache, solve_parameterized_operating_point)
from .hydraulic_resistance import MU_PA_S, ROI_TARGET_Q_M3_S
from .parameter_config import code_provenance

PARAMETERS=('s_O1','s_O2')


@dataclass(frozen=True)
class DesignEnvelope:
    lower: tuple = (.5,.5)
    upper: tuple = (2.,2.)

    def __post_init__(self):
        lo,hi=np.asarray(self.lower,float),np.asarray(self.upper,float)
        if (lo.shape!=(2,) or hi.shape!=(2,) or not np.isfinite([lo,hi]).all()
                or np.any(lo<=0) or np.any(hi<=lo)):
            raise ValueError('Explicit finite positive ordered design bounds required')

    def check(self,scales):
        s=np.asarray(scales,float)
        if (s.shape!=(2,) or not np.isfinite(s).all()
                or np.any(s<self.lower) or np.any(s>self.upper)):
            raise ValueError('Design variables outside declared finite envelope')
        return s

    def active(self,scales):
        s=self.check(scales)
        return {p:('LOWER' if abs(np.log(s[i]/self.lower[i]))<=1e-8 else
                   'UPPER' if abs(np.log(s[i]/self.upper[i]))<=1e-8 else 'NONE')
                for i,p in enumerate(PARAMETERS)}

    def as_dict(self):
        return {p:[float(self.lower[i]),float(self.upper[i])] for i,p in enumerate(PARAMETERS)}


def metrics(fractions):
    f=np.asarray(fractions,float)
    if f.shape!=(3,) or not np.isfinite(f).all():
        raise ValueError('Three finite signed fractions required')
    return dict(J_balance=float(np.sum((f-1/3)**2)),range=float(np.ptp(f)),
                std=float(np.std(f)),max_abs_deviation=float(np.max(abs(f-1/3))),
                mean_fraction=float(f.mean()))


def state_checks(state):
    positive=lambda a:bool(np.isfinite(a).all() and np.all(a>0))
    checks=dict(mass_balance=state.mass_audit['status']=='PASS',
        no_outlet_backflow=not any(state.reverse_flow_audit['port_backflow'].values()),
        finite_positive_radii=positive(state.effective_radius_m),
        finite_positive_edge_resistance=positive(state.edge_resistance_pa_s_m3))
    if not all(checks.values()):raise ValueError('Design network physical gate failed: '+str(checks))
    return checks


class DesignEvaluator:
    def __init__(self,cache,envelope):
        self.cache,self.envelope=cache,envelope
        self.trace=[]
        self.phase='baseline'

    def evaluate(self,scales):
        row=dict(evaluation=len(self.trace)+1,phase=self.phase,status='FAILED',error='')
        try:
            s=self.envelope.check(scales)
            state=solve_parameterized_operating_point(self.cache,
                ParameterizedHydraulicSpec(s_O1=float(s[0]),s_O2=float(s[1])))
            state_checks(state)
            row.update(status='PASS',s_O1=float(s[0]),s_O2=float(s[1]),
                       **metrics(state.signed_outlet_fractions),
                       mass_residual=state.mass_audit['max_relative_residual'])
            row.update({'bound_status_'+p:self.envelope.active(s)[p] for p in PARAMETERS})
            for p,q,pressure in zip(PORT_ORDER,state.port_flow_m3_s,state.port_pressure_pa):
                row['Q_'+p+'_m3_s']=float(q);row['P_'+p+'_pa']=float(pressure)
            row.update({'f_'+p:float(f) for p,f in zip(OUTLET_ORDER,state.signed_outlet_fractions)})
            self.trace.append(row)
            return state
        except Exception as exc:
            row['error']=type(exc).__name__+': '+str(exc);self.trace.append(row)
            raise

    def local(self,initial,label):
        self.phase=label
        start=len(self.trace)
        x0=np.log(self.envelope.check(initial))
        # Native bounds, never post-hoc parameter clipping.
        opt=least_squares(lambda x:self.evaluate(np.exp(x)).signed_outlet_fractions-1/3,
            x0,jac='3-point',bounds=(np.log(self.envelope.lower),np.log(self.envelope.upper)),
            ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=1000)
        state=self.evaluate(np.exp(opt.x))
        result=dict(initial_scales=list(map(float,initial)),parameters=dict(zip(PARAMETERS,np.exp(opt.x).tolist())),
            metrics=metrics(state.signed_outlet_fractions),prediction=state.summary(),
            optimizer=dict(success=bool(opt.success),status=int(opt.status),message=opt.message,
                nfev=int(opt.nfev),njev=int(opt.njev),optimality=float(opt.optimality),
                native_active_mask=opt.active_mask.tolist()),
            active_bounds=self.envelope.active(np.exp(opt.x)),physical_checks=state_checks(state),
            forward_evaluations=len(self.trace)-start)
        return result,state


def _trace_payload(rows):
    return [{k:v for k,v in row.items() if k not in ('evaluation','phase')} for row in rows]


def design_crosscheck(evaluator,grid_size=25):
    if not isinstance(grid_size,int) or grid_size<3:raise ValueError('At least 3 grid points per axis')
    baseline=evaluator.evaluate((1.,1.))
    start=len(evaluator.trace)
    first,state=evaluator.local((1.,1.),'primary_local_from_baseline')
    first_rows=_trace_payload(evaluator.trace[start:])
    evaluator.phase='coarse_log_grid';grid=[]
    for x in np.linspace(np.log(evaluator.envelope.lower[0]),np.log(evaluator.envelope.upper[0]),grid_size):
        for y in np.linspace(np.log(evaluator.envelope.lower[1]),np.log(evaluator.envelope.upper[1]),grid_size):
            evaluator.evaluate(np.exp([x,y]));grid.append(evaluator.trace[-1].copy())
    best=min(grid,key=lambda r:r['J_balance'])
    second,_=evaluator.local((best['s_O1'],best['s_O2']),'primary_local_from_grid')
    repeat_start=len(evaluator.trace)
    repeat,_=evaluator.local((1.,1.),'primary_deterministic_repeat')
    repeat_rows=_trace_payload(evaluator.trace[repeat_start:])
    pa=np.array(list(first['parameters'].values()));pb=np.array(list(second['parameters'].values()))
    cross=dict(status='CONSISTENT' if np.allclose(pa,pb,rtol=1e-6,atol=1e-8)
               and abs(first['metrics']['J_balance']-second['metrics']['J_balance'])<1e-10 else 'INCONSISTENT',
        grid_size_per_axis=grid_size,grid_points=len(grid),grid_best=best,
        parameter_absolute_difference=abs(pa-pb).tolist(),
        objective_absolute_difference=abs(first['metrics']['J_balance']-second['metrics']['J_balance']),
        grid_objective_not_better=bool(first['metrics']['J_balance']<=best['J_balance']+1e-12),
        scope='Deterministic coarse-grid and two-start crosscheck; not a mathematical global optimality proof')
    reproducible=first_rows==repeat_rows
    return baseline,first,state,dict(second_start=second,repeat=repeat,crosscheck=cross,
        deterministic_reproducibility=dict(status='PASS' if reproducible else 'FAIL',
            trace_payload_exact_equal=reproducible,compared_forward_evaluations=len(first_rows))),grid


def design_code_hashes():
    root=Path(__file__).resolve().parents[1]
    return dict(code_provenance(),balance_design=sha256(__file__),
                optimize_balanced_0d_design=sha256(root/'scripts/optimize_balanced_0d_design.py'))


def acceptance(summary):
    opt=summary['optimized'];base=summary['baseline'];cross=summary['crosscheck']
    bounds=summary['design_bounds']
    envelope=DesignEnvelope(tuple(bounds[p][0] for p in PARAMETERS),tuple(bounds[p][1] for p in PARAMETERS))
    envelope.check([opt['parameters'][p] for p in PARAMETERS])
    flags=dict(optimizer_termination=opt['optimizer']['success'],
        physical_checks=all(opt['physical_checks'].values()),
        objective_improved=opt['metrics']['J_balance']<base['metrics']['J_balance'],
        range_improved=opt['metrics']['range']<base['metrics']['range'],
        crosscheck=cross['status']=='CONSISTENT' and cross['grid_objective_not_better'],
        deterministic_reproducibility=summary['deterministic_reproducibility']['status']=='PASS',
        independent_second_start=summary['local_crosschecks']['second_start']['optimizer']['success'],
        independent_repeat=summary['local_crosschecks']['repeat']['optimizer']['success'])
    if not all(flags.values()):
        reason=('NO_BALANCE_IMPROVEMENT_WITHIN_DESIGN_ENVELOPE' if not flags['objective_improved'] or not flags['range_improved'] else 'DESIGN_VALIDATION_FAILED')
        raise ValueError(reason+': '+str(flags))
    return flags


def _seal(payload):
    return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def load_frozen_design(path):
    d=yaml.safe_load(Path(path).read_text())
    seal=d.pop('artifact_sha256',None)
    if seal!=_seal(d):raise ValueError('Frozen design integrity failure')
    if (d['status']!='FROZEN_DESIGN' or d['model_name']!=MODEL_NAME
            or d['workflow_kind']!='DESIGN_OPTIMIZATION' or d['data_identification_claim']!='NONE'
            or d['bounds_source']!='SYNTHETIC_DESIGN_ENVELOPE' or d['CFD_feedback_used'] is not False):
        raise ValueError('Invalid frozen design semantics')
    spec=ParameterizedHydraulicSpec(**d['optimized']['prediction']['parameters'])
    if (spec.mu_pa_s!=MU_PA_S or spec.distal_reference_pa!=0 or spec.o3_mode!='legacy_reference'
            or spec.terminal_resistance_O3_pa_s_m3 is not None or d['roi_target_flow_m3_s']!=ROI_TARGET_Q_M3_S):
        raise ValueError('Prescribed physical model changed')
    if d['design_bounds']!=DesignEnvelope().as_dict():raise ValueError('Only primary envelope may be frozen')
    acceptance(d)
    d['artifact_sha256']=seal
    return d


def freeze_design(summary,cache,output):
    output=Path(output)
    if output.exists():raise FileExistsError('Frozen design cannot be overwritten')
    acceptance(summary)
    if summary['provenance']['code_sha256']!=design_code_hashes():raise ValueError('Design code changed before freeze')
    if sha256(summary['provenance']['source_config_path'])!=summary['provenance']['source_config_sha256']:
        raise ValueError('Design config changed before freeze')
    if summary['input_hashes']!=dict(cache.input_hashes):raise ValueError('Geometry changed before freeze')
    if summary['design_bounds']!=DesignEnvelope().as_dict():raise ValueError('Primary envelope required')
    state=solve_parameterized_operating_point(cache,ParameterizedHydraulicSpec(**summary['optimized']['parameters']))
    state_checks(state)
    old=summary['optimized']['prediction']
    np.testing.assert_allclose(state.port_pressure_pa,[old['real_cut_pressure_pa'][p] for p in PORT_ORDER],rtol=1e-12,atol=1e-9)
    np.testing.assert_allclose(state.port_flow_m3_s,[old['port_flow_m3_s'][p] for p in PORT_ORDER],rtol=1e-12,atol=1e-28)
    data=json.loads(json.dumps(dict(summary,status='FROZEN_DESIGN',
        freeze_verification='one unchanged pure 0D forward; no optimization'),allow_nan=False))
    data['artifact_sha256']=_seal(data)
    encoded=yaml.safe_dump(data,sort_keys=False,allow_unicode=True)
    fd=os.open(output,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o444)
    with os.fdopen(fd,'w') as f:f.write(encoded);f.flush();os.fsync(f.fileno())
    return load_frozen_design(output)


def run_design(config_path,output):
    config_path=Path(config_path).resolve();cfg=yaml.safe_load(config_path.read_text())
    if (cfg['schema_version']!=1 or cfg['workflow_kind']!='DESIGN_OPTIMIZATION'
            or cfg['parameter_interpretation']!='DESIGN_VARIABLES' or cfg['data_identification_claim']!='NONE'
            or cfg['bounds_source']!='SYNTHETIC_DESIGN_ENVELOPE' or cfg['design_bounds']!=DesignEnvelope().as_dict()
            or cfg['initial_scales_dimensionless']!=[1.,1.]):raise ValueError('Invalid primary design contract')
    if (cfg['dynamic_viscosity_pa_s']!=MU_PA_S or cfg['distal_reference_pa']!=0
            or cfg['o3_mode']!='legacy_reference' or cfg['terminal_resistance_O3_pa_s_m3'] is not None
            or cfg['roi_target_flow_m3_s']!=ROI_TARGET_Q_M3_S):raise ValueError('Fixed model cannot change')
    resolve=lambda p:(config_path.parent/p).resolve()
    cache=load_geometry_cache(resolve(cfg['geometry']['graph_si_npz']),resolve(cfg['geometry']['ports_json']))
    certificate=resolve(cfg['infeasibility_certificate']);cert=json.loads(certificate.read_text())
    if cert['source_geometry_sha256']!=cache.geometry_sha256 or cert['status']!='EQUAL_SPLIT_IMPOSSIBLE_UNDER_PRESCRIBED_MODEL':
        raise ValueError('Existing infeasibility certificate does not match')
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    (output/'source_config.yaml').write_bytes(config_path.read_bytes())
    ev=DesignEvaluator(cache,DesignEnvelope())
    try:baseline,opt,state,checks,grid=design_crosscheck(ev,cfg['grid_size_per_axis'])
    finally:
        if ev.trace:write_csv(output/'design_trace.csv',ev.trace,fields=list(dict.fromkeys(k for r in ev.trace for k in r)))
    sensitivity=[]
    for upper in cfg['bound_sensitivity_upper_scales_dimensionless']:
        if upper==2.:
            run=opt
        else:
            sub=DesignEvaluator(cache,DesignEnvelope((.5,.5),(upper,upper)))
            try:run,_=sub.local((1.,1.),'bound_sensitivity')
            finally:
                if sub.trace:write_csv(output/f'bound_sensitivity_smax_{upper:g}_trace.csv',sub.trace)
        if not run['optimizer']['success']:raise ValueError('Bound sensitivity optimizer failed')
        sensitivity.append(dict(s_max_dimensionless=upper,**run['parameters'],**run['metrics'],
            **{'f_'+p:f for p,f in run['prediction']['outlet_flow_fraction'].items()},
            **{'active_'+p:v for p,v in run['active_bounds'].items()},role='0D_DIAGNOSTIC_ONLY' if upper!=2 else 'PRIMARY_DESIGN'))
    bm=metrics(baseline.signed_outlet_fractions);om=opt['metrics']
    summary=dict(schema_version=1,status='BEST_FEASIBLE_BALANCE_DESIGN_PASS',model_name=MODEL_NAME,
        workflow_kind='DESIGN_OPTIMIZATION',purpose='BEST_FEASIBLE_THREE_OUTLET_BALANCE',
        parameter_interpretation='DESIGN_VARIABLES',data_identification_claim='NONE',
        design_bounds=DesignEnvelope().as_dict(),bounds_source=cfg['bounds_source'],bounds_interpretation=cfg['bounds_interpretation'],
        equal_split_feasible=False,equal_split_required=False,infeasibility_certificate_sha256=sha256(certificate),
        ratio_bound_Q_O1_over_Q_O3=cert['ratio_bound_Q_O1_over_Q_O3'],
        theoretical_necessary_region_projection=cert['equal_sigma_least_squares_boundary_projection'],
        projection_interpretation='Larger necessary region boundary projection; not an achievable finite-envelope design',
        dynamic_viscosity_pa_s=MU_PA_S,roi_target_flow_m3_s=ROI_TARGET_Q_M3_S,
        baseline=dict(parameters=dict(s_O1=1.,s_O2=1.),metrics=bm,prediction=baseline.summary()),optimized=opt,
        optimum_classification='BOUND_ACTIVE_OPTIMUM' if any(v!='NONE' for v in opt['active_bounds'].values()) else 'INTERIOR_OPTIMUM',
        improvement=dict(delta_J=bm['J_balance']-om['J_balance'],relative_J_reduction=1-om['J_balance']/bm['J_balance'],
                         range_reduction=bm['range']-om['range'],std_reduction=bm['std']-om['std']),
        crosscheck=checks['crosscheck'],deterministic_reproducibility=checks['deterministic_reproducibility'],local_crosschecks=checks,
        bound_sensitivity=sensitivity,primary_forward_evaluations=len(ev.trace),
        input_hashes=dict(cache.input_hashes),source_geometry_sha256=cache.geometry_sha256,
        provenance=dict(source_config_sha256=sha256(config_path),source_config_path=str(config_path),
            code_sha256=design_code_hashes(),numpy_version=np.__version__,scipy_version=scipy.__version__),
        CFD_feedback_used=False,CFD_calls=0,particle_based_fitting=False,particle_RBC_calls=0)
    summary['acceptance_checks']=acceptance(summary)
    write_json(output/'design_summary.json',summary)
    write_csv(output/'coarse_grid.csv',grid)
    write_csv(output/'bound_sensitivity.csv',sensitivity)
    np.savez_compressed(output/'network_state_si.npz',pressure_pa=state.pressure_pa,edge_flow_m3_s=state.edge_flow_m3_s,
        edge_nodes=state.edge_nodes,original_node_ids=cache.node_ids,effective_radius_m=state.effective_radius_m,
        edge_resistance_pa_s_m3=state.edge_resistance_pa_s_m3,port_indices=cache.port_indices)
    return summary,cache
