"""Explicit-unit configuration and auditable command-line workflows; pure 0D."""
from dataclasses import asdict, replace
from pathlib import Path
import json
import sys
import numpy as np
import scipy
import yaml
from .audit import sha256, write_json
from .parameterized_hydraulics import (
    MODEL_NAME, PORT_ORDER, OUTLET_ORDER, ParameterizedHydraulicSpec,
    load_geometry_cache, solve_parameterized_operating_point)
from .parameter_fit import Observation, Targets, LogPrior, fit_parameters, UnidentifiableParameterizationError


def keys(mapping, allowed, required, where):
    if not isinstance(mapping,dict) or set(mapping)-set(allowed) or set(required)-set(mapping):
        raise ValueError(f'{where}: missing/unknown configuration fields; allowed={allowed}')


def number(value,where,positive=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not np.isfinite(value) or (positive and value<=0):
        raise ValueError(where+' requires an explicit finite SI number'+(' >0' if positive else ''))
    return float(value)


def code_provenance():
    names=('parameterized_hydraulics.py','parameter_fit.py','parameter_config.py','parameter_freeze.py',
           'hydraulic_resistance.py','network_solver.py','idealized_h0.py','boundary_conditions.py','audit.py')
    hashes={name:sha256(Path(__file__).parent/name) for name in names if (Path(__file__).parent/name).exists()}
    for name in ('fit_a_parameterized_0d.py','freeze_a_parameterized_0d.py'):
        path=Path(__file__).resolve().parents[1]/'scripts'/name
        hashes['scripts/'+name]=sha256(path)
    return hashes


def read_configuration(path):
    path=Path(path).resolve();cfg=yaml.safe_load(path.read_text())
    allowed=('schema_version','model','geometry','fluid','operating_point','pressure_reference','outlets','targets','optimizer','priors','regression_reference_json')
    keys(cfg,allowed,allowed[:7],'root')
    if cfg['schema_version']!=1:raise ValueError('Unsupported schema_version')
    keys(cfg['model'],('name','purpose'),('name','purpose'),'model')
    if cfg['model']['name']!=MODEL_NAME or cfg['model']['purpose'] not in ('LEGACY_REGRESSION','SYNTHETIC_VERIFICATION','EXPERIMENTAL_CALIBRATION'):
        raise ValueError('Unsupported model name/purpose')
    keys(cfg['geometry'],('graph_si_npz','ports_json'),('graph_si_npz','ports_json'),'geometry')
    keys(cfg['fluid'],('dynamic_viscosity_pa_s','density_kg_m3'),('dynamic_viscosity_pa_s','density_kg_m3'),'fluid')
    mu=number(cfg['fluid']['dynamic_viscosity_pa_s'],'mu',True)
    number(cfg['fluid']['density_kg_m3'],'density',True) # recorded, not in the steady resistance equations
    keys(cfg['operating_point'],('roi_target_flow_m3_s',),('roi_target_flow_m3_s',),'operating_point')
    q=number(cfg['operating_point']['roi_target_flow_m3_s'],'target ROI Q',True)
    keys(cfg['pressure_reference'],('distal_reference_pa','meaning'),('distal_reference_pa','meaning'),'pressure_reference')
    if cfg['pressure_reference']['meaning']!='GAUGE_REFERENCE':raise ValueError('Reference must explicitly be GAUGE_REFERENCE')
    pd=number(cfg['pressure_reference']['distal_reference_pa'],'distal_reference_pa')
    keys(cfg['outlets'],OUTLET_ORDER,OUTLET_ORDER,'outlets')
    scales={};free=[];bounds={}
    for name in ('O1','O2'):
        o=cfg['outlets'][name]
        keys(o,('downstream_model','radius_scale'),('downstream_model','radius_scale'),name)
        if o['downstream_model']!='EXISTING_GEOMETRY':raise ValueError('O1/O2 require existing downstream geometry')
        param=o['radius_scale'];keys(param,('value_dimensionless','free','parameterization','bounds_dimensionless'),('value_dimensionless','free','parameterization'),name+' radius_scale')
        if param['parameterization']!='LOG' or type(param['free']) is not bool:raise ValueError('LOG parameterization and boolean free required')
        scales['s_'+name]=number(param['value_dimensionless'],name+' radius scale',True)
        if param['free']:free.append('s_'+name)
        if 'bounds_dimensionless' in param:bounds['s_'+name]=param['bounds_dimensionless']
    o3=cfg['outlets']['O3'];keys(o3,('downstream_model','distal_pressure_pa','terminal_resistance'),('downstream_model','distal_pressure_pa'),'O3')
    if number(o3['distal_pressure_pa'],'O3 distal pressure')!=pd:
        raise ValueError('V1 uses one common distal gauge reference; heterogeneous distal pressures are not implemented')
    mode={'LEGACY_REFERENCE':'legacy_reference','TERMINAL_RESISTANCE':'terminal_resistance'}.get(o3['downstream_model'])
    if mode is None:raise ValueError('Unknown O3 model')
    resistance=None
    if mode=='terminal_resistance':
        param=o3.get('terminal_resistance')
        keys(param,('value_pa_s_m3','free','parameterization','bounds_pa_s_m3'),('value_pa_s_m3','free','parameterization'),'O3 terminal_resistance')
        if param['parameterization']!='LOG' or type(param['free']) is not bool:raise ValueError('LOG parameterization and boolean free required')
        resistance=number(param['value_pa_s_m3'],'O3 terminal resistance',True)
        if param['free']:free.append('terminal_resistance_O3_pa_s_m3')
        if 'bounds_pa_s_m3' in param:bounds['terminal_resistance_O3_pa_s_m3']=param['bounds_pa_s_m3']
    elif 'terminal_resistance' in o3:raise ValueError('Legacy O3 cannot have a resistance')
    spec=ParameterizedHydraulicSpec(**scales,mu_pa_s=mu,distal_reference_pa=pd,o3_mode=mode,terminal_resistance_O3_pa_s_m3=resistance)
    spec.validate()
    if any(p not in free for p in bounds):raise ValueError('Bounds specified for a fixed parameter')
    if cfg['model']['purpose']=='LEGACY_REGRESSION':
        if free or cfg.get('priors') or cfg.get('targets') or 'regression_reference_json' not in cfg:
            raise ValueError('Legacy regression requires fixed physics, a reference JSON and no fit targets/priors')
    opt=cfg.get('optimizer',{'method':'scipy_least_squares','max_nfev':1000})
    keys(opt,('method','max_nfev','R_ref_pa_s_m3'),('method','max_nfev'),'optimizer')
    if opt['method']!='scipy_least_squares':raise ValueError('Only scipy_least_squares supported')
    priors=[]
    for row in cfg.get('priors',[]):
        keys(row,('parameter','mean_log_dimensionless','sigma_log_dimensionless'),('parameter','mean_log_dimensionless','sigma_log_dimensionless'),'prior')
        priors.append(LogPrior(row['parameter'],number(row['mean_log_dimensionless'],'prior mean'),number(row['sigma_log_dimensionless'],'prior sigma',True)))
    resolve=lambda p:(path.parent/p).resolve()
    return dict(config=cfg,path=path,config_sha256=sha256(path),spec=spec,target_q=q,free=tuple(free),
                bounds=bounds,priors=tuple(priors),optimizer=opt,
                geometry_path=resolve(cfg['geometry']['graph_si_npz']),ports_path=resolve(cfg['geometry']['ports_json']))


def configuration_targets(loaded,cache):
    cfg=loaded['config'];t=cfg.get('targets',{})
    allowed=('kind','evidence_reference','synthetic_truth','sigma_flow_fraction_dimensionless',
             'outlet_flow_fraction','mean_outlet_flow_m3_s','real_cut_pressure_pa','pressure_difference_pa')
    keys(t,allowed,('kind','evidence_reference'),'targets')
    if t['kind'] not in ('SYNTHETIC_TARGET','EXPERIMENTAL_TARGET'):raise ValueError('CFD/particle target sources forbidden')
    purpose=cfg['model']['purpose']
    if (purpose=='SYNTHETIC_VERIFICATION')!=(t['kind']=='SYNTHETIC_TARGET') or purpose=='LEGACY_REGRESSION':
        raise ValueError('Model purpose and target kind disagree')
    obs=[];truth=None
    if 'synthetic_truth' in t:
        if t['kind']!='SYNTHETIC_TARGET':raise ValueError('Synthetic truth cannot be experimental data')
        if any(k in t for k in allowed[4:]):raise ValueError('Choose synthetic generation or explicit targets, not both')
        truth=t['synthetic_truth'];keys(truth,('s_O1_dimensionless','s_O2_dimensionless','terminal_resistance_O3_pa_s_m3'),('s_O1_dimensionless','s_O2_dimensionless'),'synthetic_truth')
        kw={k:number(v,k,True) for k,v in truth.items()}
        kw['s_O1']=kw.pop('s_O1_dimensionless');kw['s_O2']=kw.pop('s_O2_dimensionless')
        true_spec=replace(loaded['spec'],**kw)
        r=solve_parameterized_operating_point(cache,true_spec,loaded['target_q'])
        sigma=number(t['sigma_flow_fraction_dimensionless'],'fraction sigma',True)
        obs=[Observation('flow_fraction',p,float(f),sigma) for p,f in zip(OUTLET_ORDER,r.signed_outlet_fractions)]
        truth=dict(parameters=asdict(true_spec),prediction=r.summary(),source='SYNTHETIC_PURE_0D_FORWARD')
    else:
        if 'sigma_flow_fraction_dimensionless' in t:raise ValueError('Explicit observations each require their own sigma')
        for key,kind,unit in [('outlet_flow_fraction','flow_fraction','dimensionless'),
                              ('mean_outlet_flow_m3_s','outlet_flow_m3_s','m3_s'),
                              ('real_cut_pressure_pa','real_cut_pressure_pa','pa')]:
            for port,row in t.get(key,{}).items():
                fields=('value_'+unit,'sigma_'+unit);keys(row,fields,fields,key)
                obs.append(Observation(kind,port,number(row[fields[0]],fields[0]),number(row[fields[1]],fields[1],True)))
        for row in t.get('pressure_difference_pa',[]):
            fields=('from_port','to_port','value_pa','sigma_pa');keys(row,fields,fields,'pressure_difference_pa')
            obs.append(Observation('pressure_difference_pa',row['from_port'],number(row['value_pa'],'pressure difference'),number(row['sigma_pa'],'pressure sigma',True),row['to_port']))
    target=Targets(t['kind'],tuple(obs),t['evidence_reference']).canonical()
    return target,truth


def run_configuration(config_path,output):
    loaded=read_configuration(config_path)
    cache=load_geometry_cache(loaded['geometry_path'],loaded['ports_path'])
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    (output/'source_config.yaml').write_bytes(loaded['path'].read_bytes())
    provenance=dict(source_config_sha256=loaded['config_sha256'],source_config_path=str(loaded['path']),
                    source_geometry_path=str(loaded['geometry_path']),source_ports_path=str(loaded['ports_path']),
                    code_sha256=code_provenance(),numpy_version=np.__version__,scipy_version=scipy.__version__,
                    pyyaml_version=yaml.__version__,python_version=sys.version)
    if loaded['free']:
        targets,truth=configuration_targets(loaded,cache)
        if truth is not None:write_json(output/'synthetic_truth.json',truth)
        try:
            result=fit_parameters(cache,loaded['spec'],targets,free_parameters=loaded['free'],
                target_roi_flow_m3_s=loaded['target_q'],priors=loaded['priors'],physical_bounds=loaded['bounds'],
                max_nfev=loaded['optimizer']['max_nfev'],R_ref_pa_s_m3=loaded['optimizer'].get('R_ref_pa_s_m3',1e17),
                trace_path=output/'parameter_trace.csv',provenance=provenance,state_output_path=output/'network_state_si.npz')
        except UnidentifiableParameterizationError as exc:
            write_json(output/'fit_rejection.json',dict(status='REJECTED_UNIDENTIFIABLE',audit=exc.audit,provenance=provenance,CFD_calls=0,particle_simulation_calls=0))
            raise
    else:
        if loaded['config']['model']['purpose']!='LEGACY_REGRESSION':
            raise ValueError('No free parameters: use explicit LEGACY_REGRESSION purpose in this CLI')
        state=solve_parameterized_operating_point(cache,loaded['spec'],loaded['target_q'])
        reference=(loaded['path'].parent/loaded['config']['regression_reference_json']).resolve()
        reference_source='CURRENT_SAVED_H0_0D_REPORT'
        if not reference.exists():
            reference=Path(__file__).resolve().parents[1]/'configs/parameterized_hydraulics/legacy_h0_snapshot.json'
            reference_source='RECORDED_CURRENT_SNAPSHOT_FALLBACK_NOT_EXPERIMENTAL_TRUTH'
        gold=json.loads(reference.read_text()) # 0D-only H0 regression oracle; never a fit target
        pressures=np.array([r['pressure_realcut_Pa'] for r in gold['ports']])
        fractions=np.array([r['signed_fraction'] for r in gold['ports'][1:]])
        np.testing.assert_allclose(state.port_pressure_pa,pressures,rtol=2e-11,atol=1e-9)
        np.testing.assert_allclose(state.signed_outlet_fractions,fractions,rtol=2e-11,atol=2e-13)
        np.testing.assert_allclose(state.port_flow_m3_s[0],gold['ports'][0]['signed_Q_m3s'],rtol=2e-14,atol=0)
        if loaded['spec']!=ParameterizedHydraulicSpec():raise ValueError('Legacy regression config must use exact default H0 physics')
        result=dict(status='FORWARD_VALIDATED',model_name=MODEL_NAME,parameter_source='LEGACY_REFERENCE_REGRESSION_NOT_EXPERIMENT',
            initial_parameters=asdict(loaded['spec']),final_parameters=asdict(loaded['spec']),
            roi_target_flow_m3_s=loaded['target_q'],prediction=state.summary(),provenance=provenance,
            source_geometry_sha256=cache.geometry_sha256,input_hashes=dict(cache.input_hashes),
            regression=dict(reference_source=reference_source,reference_sha256=sha256(reference),maximum_pressure_error_pa=float(np.max(abs(state.port_pressure_pa-pressures))),
                            maximum_fraction_error_dimensionless=float(np.max(abs(state.signed_outlet_fractions-fractions))),
                            inlet_flow_error_m3_s=float(state.port_flow_m3_s[0]-gold['ports'][0]['signed_Q_m3s']),
                            legacy_mass_audit=gold['mass']),CFD_calls=0,particle_simulation_calls=0)
    write_json(output/'fit_summary.json',result)
    return result
