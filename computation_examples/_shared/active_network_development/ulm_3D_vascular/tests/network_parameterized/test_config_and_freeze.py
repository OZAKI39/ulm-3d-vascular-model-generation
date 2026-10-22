from pathlib import Path
import json
import os
import yaml
import pytest
from network_1d0d.parameter_config import read_configuration,run_configuration
from network_1d0d.parameter_freeze import freeze_parameters,load_frozen_artifact

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/'configs/parameterized_hydraulics/synthetic_radius_fit_example.yaml'


def local_config(tmp_path):
    c=yaml.safe_load(CONFIG.read_text())
    for k,v in c['geometry'].items():c['geometry'][k]=str((CONFIG.parent/v).resolve())
    p=tmp_path/'config.yaml';p.write_text(yaml.safe_dump(c));return p,c


def test_freeze_artifact_roundtrip(tmp_path):
    config,_=local_config(tmp_path)
    result=run_configuration(config,tmp_path/'fit')
    output=tmp_path/'frozen_parameterized_boundary.yaml'
    frozen=freeze_parameters(config,tmp_path/'fit/fit_summary.json',output)
    assert frozen==load_frozen_artifact(output)
    assert frozen['parameters']==result['final_parameters']
    assert frozen['status']=='FROZEN' and frozen['feedback_from_CFD'] is False
    assert output.stat().st_mode & 0o222 == 0
    assert frozen['0d_prediction']['outlet_flow_fraction']==result['prediction']['outlet_flow_fraction']
    with pytest.raises(FileExistsError):freeze_parameters(config,tmp_path/'fit/fit_summary.json',output)
    tampered=yaml.safe_load(output.read_text());tampered['parameters']['s_O1']=1.1
    other=tmp_path/'tampered.yaml';other.write_text(yaml.safe_dump(tampered))
    with pytest.raises(ValueError,match='integrity'):load_frozen_artifact(other)


def test_changed_config_cannot_be_frozen(tmp_path):
    config,c=local_config(tmp_path);run_configuration(config,tmp_path/'fit')
    c['fluid']['dynamic_viscosity_pa_s']*=2;config.write_text(yaml.safe_dump(c))
    with pytest.raises(ValueError,match='config changed'):freeze_parameters(config,tmp_path/'fit/fit_summary.json',tmp_path/'frozen.yaml')


def test_failed_fit_cannot_be_frozen(tmp_path):
    config,_=local_config(tmp_path);run_configuration(config,tmp_path/'fit')
    file=tmp_path/'fit/fit_summary.json';summary=json.loads(file.read_text());summary['status']='FIT_NOT_CONVERGED';file.write_text(json.dumps(summary))
    with pytest.raises(ValueError,match='converged'):freeze_parameters(config,file,tmp_path/'frozen.yaml')


def test_workflow_never_overwrites_prior_run(tmp_path):
    config,_=local_config(tmp_path);run_configuration(config,tmp_path/'fit')
    with pytest.raises(FileExistsError):run_configuration(config,tmp_path/'fit')


def test_unknown_unit_field_is_rejected(tmp_path):
    config,c=local_config(tmp_path);c['fluid']['dynamic_viscosity']=c['fluid'].pop('dynamic_viscosity_pa_s');config.write_text(yaml.safe_dump(c))
    with pytest.raises(ValueError,match='fields'):read_configuration(config)


def test_heterogeneous_distal_pressures_are_not_silently_scaled(tmp_path):
    config,c=local_config(tmp_path);c['outlets']['O3']['distal_pressure_pa']=50;config.write_text(yaml.safe_dump(c))
    with pytest.raises(ValueError,match='common distal'):read_configuration(config)


def test_legacy_configuration_runs_without_optimizer(tmp_path,monkeypatch):
    import network_1d0d.parameter_config as pc
    def forbidden(*a,**k):raise AssertionError('Legacy regression must not fit')
    monkeypatch.setattr(pc,'fit_parameters',forbidden)
    r=run_configuration(ROOT/'configs/parameterized_hydraulics/legacy_h0_regression.yaml',tmp_path/'legacy')
    assert r['status']=='FORWARD_VALIDATED'
    assert r['regression']['maximum_pressure_error_pa']<1e-9


def test_documented_legacy_snapshot_fallback(tmp_path):
    original=ROOT/'configs/parameterized_hydraulics/legacy_h0_regression.yaml'
    c=yaml.safe_load(original.read_text())
    for k,v in c['geometry'].items():c['geometry'][k]=str((original.parent/v).resolve())
    c['regression_reference_json']='missing_historical_report.json'
    path=tmp_path/'config.yaml';path.write_text(yaml.safe_dump(c))
    r=run_configuration(path,tmp_path/'fallback')
    assert r['status']=='FORWARD_VALIDATED'
    assert r['regression']['reference_source'].startswith('RECORDED_CURRENT_SNAPSHOT_FALLBACK')
