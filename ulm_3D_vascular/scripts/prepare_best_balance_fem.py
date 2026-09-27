"""Frozen primary design -> regressed one-way handoff -> ONE fixed FEM case.

This command creates inputs only; it never launches any solver.
"""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from network_1d0d.audit import sha256,write_json
from network_1d0d.balance_design import load_frozen_design
from network_1d0d.parameterized_fem_handoff import export_frozen_design,transfer_prediction
from network_1d0d.fem_h0_case import create_case

CASE_NAME='mean-2p0-mmps-A-best-feasible-balance-v1'


def prepare(frozen_path,source,target,report):
    frozen_path,source,target,report=map(Path,(frozen_path,source,target,report))
    f=load_frozen_design(frozen_path)
    if target.name!=CASE_NAME:raise ValueError('Only the single primary scientific case is allowed')
    if target.exists():raise FileExistsError('Never overwrite/create a second case')
    report.mkdir(parents=True,exist_ok=True)
    registry=report/'scientific_case_registry.json'
    if registry.exists():raise FileExistsError('A scientific case is already registered')
    reference=ROOT/'reports/a_network_1d0d_boundary_v2_idealized/data'
    provenance=json.loads((reference/'3d_case_config_diff.json').read_text())
    if sha256(source/'run/solver.xml')!=provenance['new_solver_sha256']:
        raise ValueError('Source XML does not match authoritative H0 case provenance')
    for name,digest in provenance['unchanged_input_hashes'].items():
        if sha256(source/name)!=digest:raise ValueError('Authoritative source input changed: '+name)
    ports=ROOT/'reports/a_network_1d0d_boundary_v1/data/roi_ports_in_a.json'
    handoff=export_frozen_design(frozen_path,source,ports)
    # Regression uses the newly computed full-A baseline, not CFD flow.
    baseline=transfer_prediction(f['baseline']['prediction'],handoff['extension_measurements'])
    old=json.loads((reference/'roi_fixed_pressure_bc_H0.json').read_text())
    errors={}
    for field in ('R_extension_Pa_s_m3','pressure_cap_raw_Pa','pressure_cap_shifted_Pa'):
        aa=np.array([r[field] for r in baseline['ports']]);bb=np.array([r[field] for r in old['ports']])
        np.testing.assert_allclose(aa,bb,rtol=2e-12,atol=1e-9)
        errors[field]=float(abs(aa-bb).max())
    regression=dict(status='PASS',maximum_absolute_errors=errors,reference_sha256=sha256(reference/'roi_fixed_pressure_bc_H0.json'),
                    method='Full-A s1=s2=1 baseline + current extension geometry only; CFD calls=0')
    write_json(report/'handoff_H0_regression.json',regression)
    handoff_path=report/'best_feasible_fem_handoff.json';write_json(handoff_path,handoff)
    diff=create_case(source,target,handoff_path)
    write_json(target/'configuration_diff.json',diff)
    write_json(report/'configuration_diff.json',diff)
    inputs=json.loads((target/'input_hashes.json').read_text())
    # Pin new case, frozen design and provenance before any external process exists.
    for name,path in [('frozen_balance_design.yaml',frozen_path),('best_feasible_fem_handoff.json',handoff_path)]:
        (target/name).write_bytes(path.read_bytes())
    preflight=dict(status='PASS',starting_commit='c5635dc0b0a2404f8eba9b347784d0429071f795',
        workflow_kind='DESIGN_OPTIMIZATION',design_bounds=f['design_bounds'],best_parameters=f['optimized']['parameters'],
        baseline=f['baseline'],optimized=f['optimized'],improvement=f['improvement'],
        design_crosscheck=f['crosscheck'],deterministic_reproducibility=f['deterministic_reproducibility'],
        frozen_design_sha256=sha256(frozen_path),handoff_sha256=sha256(handoff_path),
        extension_R_pa_s_m3={r['port']:r['R_extension_Pa_s_m3'] for r in handoff['ports']},
        raw_cap_pressure_pa={r['port']:r['pressure_cap_raw_Pa'] for r in handoff['ports']},
        applied_cap_pressure_pa={r['port']:r['pressure_cap_shifted_Pa'] for r in handoff['ports']},
        input_hashes=inputs,configuration_diff=diff,handoff_regression=regression,
        source_case_provenance_sha256=sha256(reference/'3d_case_config_diff.json'),
        case_name=CASE_NAME,local_case=str(target.resolve()),
        solver_sha256='0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7',
        petsc_library_sha256='b79e98d278e105aa8a1848bd08d66ef36a3e2425458610f97213d3fc8e45c5ef',
        runner_sha256=sha256(ROOT/'scripts/solve_a_h0_fem_remote.py'),
        MPI_ranks=1,OMP_NUM_THREADS=1,CUDA_VISIBLE_DEVICES='0',
        planned_CFD_run_count=1,CFD_feedback_used=False,particle_RBC_calls=0)
    write_json(target/'final_preflight.json',preflight);write_json(report/'final_preflight.json',preflight)
    write_json(registry,dict(scientific_case_count=1,case_name=CASE_NAME,local_case=str(target.resolve()),
                            solver_launches=0,planned_CFD_run_count=1,parameter_retuning_permitted=False))
    return preflight


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('frozen-design','source-case','target-case','report'):p.add_argument('--'+arg,type=Path,required=True)
    a=p.parse_args();r=prepare(a.frozen_design,a.source_case,a.target_case,a.report)
    print(json.dumps({k:r[k] for k in ('status','design_bounds','best_parameters','extension_R_pa_s_m3','raw_cap_pressure_pa','applied_cap_pressure_pa','input_hashes','planned_CFD_run_count')},indent=2))
