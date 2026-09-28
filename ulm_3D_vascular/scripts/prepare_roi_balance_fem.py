"""Frozen ROI-only design -> current extension geometry -> one controlled case."""
from pathlib import Path
import argparse
import json
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from network_1d0d.audit import sha256,write_json,write_csv
from network_1d0d.roi_boundary_design import load_frozen_design
from network_1d0d.parameterized_fem_handoff import measure_extension_geometry,transfer_prediction
from network_1d0d.fem_h0_case import create_case,serialize_fixed_pressure_xml,audit_xml_change

CASE_NAME='mean-2p0-mmps-A-ROI-only-balanced-pressure-v1'


def prepare(frozen_path,source,target,report):
    frozen_path,source,target,report=map(Path,(frozen_path,source,target,report))
    f=load_frozen_design(frozen_path)
    if f['regression_status']!='PASS' or f['mass_audit']['status']!='PASS':raise ValueError('ROI gate failed')
    if target.name!=CASE_NAME or target.exists():raise ValueError('Only one new case, without overwriting')
    registry=report/'scientific_case_registry.json'
    if registry.exists():raise FileExistsError('Scientific case already registered')
    reference=ROOT/'reports/a_network_1d0d_boundary_v2_idealized/data'
    provenance=json.loads((reference/'3d_case_config_diff.json').read_text())
    if sha256(source/'run/solver.xml')!=provenance['new_solver_sha256']:raise ValueError('Source must be authoritative H0')
    for name,digest in provenance['unchanged_input_hashes'].items():
        if sha256(source/name)!=digest:raise ValueError('Source input changed: '+name)
    ports_path=ROOT/'reports/a_network_1d0d_boundary_v1/data/roi_ports_in_a.json'
    if sha256(ports_path)!=f['port_mapping_sha256']:raise ValueError('Port map changed')
    mapping=json.loads(ports_path.read_text())
    ext=measure_extension_geometry(source,mapping['ports'],mu_pa_s=f['dynamic_viscosity_pa_s'])
    handoff=transfer_prediction(dict(real_cut_pressure_pa=f['gauge_realcut_pressure_pa'],port_flow_m3_s=f['outlet_flow_m3_s']),ext)
    handoff.update(source_model=f['model_name'],source_frozen_design_sha256=sha256(frozen_path),
        full_A_network_used_for_boundary_generation=False,extension_measurements=ext,
        ROI_realcut_gauge_added_pa=f['gauge_shift_pa'],FEM_cap_gauge_added_pa=handoff['gauge_added_constant_Pa'],
        reference_pressure_definition='ROI-only mathematical outlet O3=0 -> common real-cut shift -> independent common cap shift',
        CFD_feedback_used=False,mu_pa_s=f['dynamic_viscosity_pa_s'])
    handoff_path=report/'roi_fem_handoff.json';write_json(handoff_path,handoff)
    write_csv(report/'roi_fem_cap_pressures.csv',handoff['ports'])
    values=[r['pressure_cap_shifted_Pa'] for r in handoff['ports']]
    xml=(source/'run/solver.xml').read_bytes()
    audit_xml_change(xml,serialize_fixed_pressure_xml(xml,values),values)
    # Source and design gates run before the sole target directory is created.
    diff=create_case(source,target,handoff_path)
    write_json(target/'configuration_diff.json',diff);write_json(report/'configuration_diff.json',diff)
    for name,path in [('roi_boundary_design_frozen.yaml',frozen_path),('roi_fem_handoff.json',handoff_path)]:
        (target/name).write_bytes(path.read_bytes())
    pre=dict(status='PASS',starting_commit='69c6af94eabb87ba6db0823134c4b92b30a069db',
        workflow_kind='ROI_BOUNDARY_DESIGN',model_name=f['model_name'],ROI_design_status=f['status'],
        full_A_network_used_for_boundary_generation=False,regression_status=f['regression_status'],
        frozen_design_sha256=sha256(frozen_path),handoff_sha256=sha256(handoff_path),
        applied_cap_pressure_pa=dict(zip(('O1','O2','O3'),values)),input_hashes=json.loads((target/'input_hashes.json').read_text()),
        configuration_diff=diff,case_name=CASE_NAME,local_case=str(target.resolve()),
        solver_sha256='0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7',
        petsc_library_sha256='b79e98d278e105aa8a1848bd08d66ef36a3e2425458610f97213d3fc8e45c5ef',
        runner_sha256=sha256(ROOT/'scripts/solve_a_h0_fem_remote.py'),
        MPI_ranks=1,OMP_NUM_THREADS=1,CUDA_VISIBLE_DEVICES='0',planned_CFD_run_count=1,CFD_feedback_used=False,particle_RBC_calls=0)
    write_json(target/'final_preflight.json',pre);write_json(report/'final_preflight.json',pre)
    with registry.open('x') as stream:json.dump(dict(scientific_case_count=1,case_name=CASE_NAME,local_case=str(target.resolve()),planned_CFD_run_count=1,parameter_retuning_permitted=False),stream,indent=2)
    return pre


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('frozen-design','source-case','target-case','report'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();r=prepare(a.frozen_design,a.source_case,a.target_case,a.report)
    print(json.dumps({k:r[k] for k in ('status','applied_cap_pressure_pa','planned_CFD_run_count')},indent=2))
