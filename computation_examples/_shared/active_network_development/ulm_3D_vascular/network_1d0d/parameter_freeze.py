"""Explicit, no-overwrite, hash-sealed parameter freeze; no optimizer or CFD."""
from dataclasses import asdict
from pathlib import Path
import hashlib
import json
import os
import numpy as np
import yaml

from .audit import sha256
from .parameterized_hydraulics import (
    MODEL_NAME, PORT_ORDER, OUTLET_ORDER, ParameterizedHydraulicSpec,
    load_geometry_cache, solve_parameterized_operating_point)
from .parameter_config import read_configuration, code_provenance


def _digest(payload):
    return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def load_frozen_artifact(path):
    data=yaml.safe_load(Path(path).read_text())
    if not isinstance(data,dict):raise ValueError('Invalid frozen artifact')
    seal=data.pop('artifact_sha256',None)
    if not isinstance(seal,str) or _digest(data)!=seal:raise ValueError('Frozen artifact integrity check failed')
    if data.get('status')!='FROZEN' or data.get('model_name')!=MODEL_NAME or data.get('schema_version')!=1:
        raise ValueError('Wrong frozen schema/model/status')
    spec=ParameterizedHydraulicSpec(**data['parameters']);spec.validate()
    if data['feedback_from_CFD'] is not False:raise ValueError('CFD feedback prohibited')
    if data['parameter_source'] not in ('SYNTHETIC_TARGET','EXPERIMENTAL_TARGET','LEGACY_REFERENCE_REGRESSION_NOT_EXPERIMENT'):
        raise ValueError('Prohibited frozen parameter source')
    data['artifact_sha256']=seal
    return data


def freeze_parameters(config_path,fit_summary_path,output_path):
    """Recompute one pure-0D prediction at the fixed optimum, never refit.

    Read-only mode plus an exclusive file creation prevents accidental
    overwrites. The canonical hash detects editing; OS owners can still alter
    files, so this is a verifiable immutable artifact contract, not a security
    guarantee against a malicious filesystem owner.
    """
    output=Path(output_path)
    if output.exists():raise FileExistsError('Frozen artifact cannot be overwritten: '+str(output))
    loaded=read_configuration(config_path)
    fit=json.loads(Path(fit_summary_path).read_text())
    if fit['status'] not in ('FIT_CONVERGED','FORWARD_VALIDATED'):
        raise ValueError('Only a validated forward or converged fit may be frozen')
    provenance=fit['provenance']
    if loaded['config_sha256']!=provenance['source_config_sha256']:
        raise ValueError('Source config changed after fitting')
    if code_provenance()!=provenance['code_sha256']:
        raise ValueError('Model code changed after fitting; produce a new audited fit')
    if fit.get('CFD_calls')!=0 or fit.get('particle_simulation_calls')!=0:
        raise ValueError('CFD/particle-assisted result cannot be frozen')
    cache=load_geometry_cache(loaded['geometry_path'],loaded['ports_path'])
    if cache.geometry_sha256!=fit['source_geometry_sha256'] or dict(cache.input_hashes)!=fit['input_hashes']:
        raise ValueError('Source geometry/mapping changed after fitting')
    if fit['roi_target_flow_m3_s']!=loaded['target_q']:
        raise ValueError('Operating point changed after fitting')
    if fit['status']=='FIT_CONVERGED':
        audit=fit['identifiability_audit']
        if not fit['optimizer_status']['success'] or audit['final_augmented_jacobian']['rank']!=len(loaded['free']):
            raise ValueError('Convergence/identifiability gate failed')
    final=ParameterizedHydraulicSpec(**fit['final_parameters'])
    # A fit is allowed to change only the explicitly free parameters.
    for name,value in asdict(loaded['spec']).items():
        if name not in loaded['free'] and getattr(final,name)!=value:
            raise ValueError('Nonfree physical parameter changed: '+name)
    state=solve_parameterized_operating_point(cache,final,loaded['target_q'])
    old=fit['prediction']
    np.testing.assert_allclose(state.port_pressure_pa,[old['real_cut_pressure_pa'][p] for p in PORT_ORDER],rtol=1e-12,atol=1e-9)
    np.testing.assert_allclose(state.port_flow_m3_s,[old['port_flow_m3_s'][p] for p in PORT_ORDER],rtol=1e-12,atol=1e-28)
    prediction=dict(inlet_pressure_pa=float(state.port_pressure_pa[0]),
        inlet_flow_m3_s=float(state.port_flow_m3_s[0]),
        outlet_pressure_pa=dict(zip(OUTLET_ORDER,state.port_pressure_pa[1:].tolist())),
        outlet_flow_m3_s=dict(zip(OUTLET_ORDER,state.port_flow_m3_s[1:].tolist())),
        outlet_flow_fraction=dict(zip(OUTLET_ORDER,state.signed_outlet_fractions.tolist())),
        mass_residual_dimensionless=state.mass_audit['max_relative_residual'])
    data=dict(schema_version=1,model_name=MODEL_NAME,status='FROZEN',parameter_source=fit['parameter_source'],
        source_config_sha256=loaded['config_sha256'],source_geometry_sha256=cache.geometry_sha256,
        source_geometry_file_sha256=cache.input_hashes['geometry_npz_sha256'],
        source_port_mapping_sha256=cache.input_hashes['port_mapping_sha256'],source_fit_summary_sha256=sha256(fit_summary_path),
        source_code_sha256=code_provenance(),dynamic_viscosity_pa_s=final.mu_pa_s,
        density_kg_m3=loaded['config']['fluid']['density_kg_m3'],roi_target_flow_m3_s=loaded['target_q'],
        pressure_reference=dict(distal_reference_pa=final.distal_reference_pa,meaning='GAUGE_REFERENCE'),
        parameters=asdict(final),feedback_from_CFD=False,particle_based_fitting=False,
        freeze_verification='one pure-0D forward at fixed final parameters; no optimization',
        **{'0d_prediction':prediction})
    for name in ('O1','O2'):
        data[name]=dict(downstream_model='EXISTING_GEOMETRY',radius_scale_dimensionless=getattr(final,'s_'+name),
                        source='EXISTING_GEOMETRY_EFFECTIVE_RADIUS_PARAMETER',parameter_source=fit['parameter_source'],real_cut_radius_scaled=False)
    data['O3']=dict(downstream_model=final.o3_mode.upper(),distal_pressure_pa=final.distal_reference_pa,
                    terminal_resistance_pa_s_m3=final.terminal_resistance_O3_pa_s_m3,
                    source='MODEL_CLOSURE_PARAMETER' if final.o3_mode=='terminal_resistance' else 'DIRICHLET_REFERENCE_TERMINAL',
                    augmented_network_metadata=state.closure_provenance)
    data['artifact_sha256']=_digest(data)
    output.parent.mkdir(parents=True,exist_ok=True)
    encoded=yaml.safe_dump(data,sort_keys=False,allow_unicode=True).encode()
    fd=os.open(output,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o444)
    with os.fdopen(fd,'wb') as f:f.write(encoded);f.flush();os.fsync(f.fileno())
    return load_frozen_artifact(output)
