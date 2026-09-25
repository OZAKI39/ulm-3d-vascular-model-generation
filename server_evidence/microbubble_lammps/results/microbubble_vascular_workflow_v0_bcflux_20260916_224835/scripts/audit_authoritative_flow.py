#!/usr/bin/env python3
"""Trace the saved physical Q, including the documented supersession of s3/s4 Q."""
from pathlib import Path
import hashlib,json,shutil,subprocess,os
import yaml
S=Path(__file__).resolve().parents[1]
V=Path('/home/lzy/projects/ulm_3D_vascular')
P=Path('/home/lzy/projects/compre_output/pure_fluid_new_medium_smoke/20260915_161937')
T=Path('/home/lzy/projects/compre_output/vast4090_stage4_restore/20260915_150938/stage4_200')
R=V/'outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628'
G=V/'outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611'
E=S/'provenance/flow_bc_evidence';E.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
records=[]
def read(p,label,parse=True):
    dest=E/(label+p.suffix);shutil.copy2(p,dest)
    records.append(dict(label=label,source_file=str(p),source_run=str(p.parent),sha256=sha(p),saved_copy=str(dest.relative_to(S))))
    if not parse:return p.read_text()
    return yaml.safe_load(p.read_text())
network=read(R/'global_1d/solver_summary.json','global_1d_solver')
roi=read(R/'roi/boundary_conditions.json','roi_boundary_conditions')
read(R/'roi/port_classification.csv','roi_port_transfer',False)
read(R/'qc/run_summary.json','roi_run_summary')
original=read(G/'bc/boundary_conditions_original.json','surface_original_bc')
refined=read(G/'bc/boundary_conditions_vmtk_boundarynormal_crossseam.json','surface_final_bc')
read(G/'bc/extension_pressure_correction_vmtk_boundarynormal_crossseam.csv','surface_transfer',False)
read(G/'boundaries/boundary_manifest.csv','surface_boundary_manifest',False)
config=read(V/'configs/cfd_flow.yaml','current_cfd_config')
accepted=read(V/'outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/qc/reference_scaled_base_runtime_contract.json','accepted_cfd_runtime')
physical=read(T/'contracts/physical_bc_contract.json','stage4_physical_bc')
lattice=read(T/'contracts/lattice_unit_contract.json','stage4_lattice')
inp=read(P/'inputs/NUMERICS_INPUT.json','pbs_bsa_input')
num=read(P/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json','pbs_bsa_contract')
params=read(P/'contracts/solver_parameters.txt','pbs_bsa_actual_parameters',False)
read(P/'inputs/OLD_SOLVER_PARAMETERS.txt','stage4_parameters',False)
read(P/'provenance/SOLVER_PARAMETER_GENERATION.json','parameter_generation')
read(P/'provenance/NUMERICS_GENERATION_RECEIPT.json','numerics_generation')
execution=read(P/'provenance/EXECUTION_SHA256.json','execution_hashes')
read(P/'scripts/prepare_numerics.py','pbs_bsa_generator',False)
read(P/'source/vascularPoC.cpp','palabos_case_source',False)
actual=read(P/'run_5000/diagnostics/ACTUAL_NUMERICS_READBACK.json','actual_readback')
cmd=read(P/'run_5000/diagnostics/inlet_command.json','actual_inlet_command')
read(P/'run_5000/RUN_STARTED.json','palabos_run_started')
read(P/'run_5000/RUN_TERMINAL.json','palabos_run_terminal')
Q=config['boundary_conditions']['target_volume_flow_m3_s'];oldQ=roi['inlet']['flow_rate_m3_s']
checks={
 'roi_Q_preserved_in_surface_original':oldQ==original['inlet']['flow_rate_m3_s'],
 'roi_Q_preserved_in_surface_final':oldQ==refined['inlet']['flow_rate_m3_s'],
 'same_port_through_surface':roi['inlet']['port_id']==refined['inlet']['port_id'],
 'current_Q_matches_accepted_CFD':Q==accepted['contract']['target_volume_flow_m3_s'],
 'current_Q_matches_stage4':Q==physical['inlet']['target_volume_flow_m3_s']==lattice['target_Q_phys'],
 'supersession_explicit':'supersedes earlier s3/s4 Q' in json.dumps(physical),
 'current_config_hash_matches_saved_evidence':any(e.get('source_sha256')==sha(V/'configs/cfd_flow.yaml') for e in physical['evidence']),
 'PBS_BSA_Q_preserved':Q==inp['Qtarget_m3_s']==actual['Qtarget']==cmd['physical_Qtarget']==float(params.splitlines()[2].split()[10]),
 'PBS_actual_contract_hash':actual['numerics_contract_sha256']==sha(P/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json'),
 'PBS_runtime_input_hashes':all(sha(P/n)==execution[n] for n in ['contracts/solver_parameters.txt','contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json','inputs/NUMERICS_INPUT.json','inputs/OLD_SOLVER_PARAMETERS.txt','source/vascularPoC.cpp','scripts/prepare_numerics.py']),
 'nominal_integral_error':abs(cmd['nominal_integral_error'])<1e-12,
 'numerical_multiplier_separate':abs(cmd['numerical_profile_integral_m3_s']/Q-cmd['multiplier'])<1e-12,
}
f=json.loads((S/'contracts/FROZEN_FLOW_CONTRACT.json').read_text())
checks['frozen_snapshot_hash']=sha(Path(f['source_file']))==f['source_sha256']
checks['frozen_field_hash']=sha(S/f['field_file'])==f['field_sha256']
checks['field_same_PBS_contract']=f['PBS_BSA_numerics_sha256']==actual['numerics_contract_sha256']
audit=dict(status='PASS' if all(checks.values()) else 'BLOCKED',checks=checks,
 source_file=str(P/'contracts/solver_parameters.txt'),source_run=str(P/'run_5000'),
 boundary_index=0,port_id=roi['inlet']['port_id'],entity_id=4,global_edge_id=814,
 authoritative_flow_rate=Q,unit='m^3/s',sign='positive into domain, -u dot outward normal',
 historical_roi_flow_m3_s=oldQ,preserved_during_CFD_surface_preparation=True,
 historical_roi_Q_used_unchanged_by_PBS_BSA=False,
 supersession='Saved current CFD target explicitly supersedes earlier s3/s4 Q; accepted downstream target -> Stage4 -> PBS/BSA preserves current Q.',
 authoritative_current_Q_used_by_PBS_BSA=True,
 profile_lineage='ROI PARABOLIC -> current CFD adaptive flux pressure -> Palabos uniform normal nominal Q/native cap area with ramp and numerical multiplier',
 physical_target_m3_s=Q,numerical_multiplier=cmd['multiplier'],numerical_profile_integral_m3_s=cmd['numerical_profile_integral_m3_s'],
 multiplier_validation='NOT_PERFORMED_FOR_NEW_MEDIUM',flow_physics_status='ENGINEERING_TRANSIENT_FIELD_ONLY',
 source_is_experimental_measurement=False,flow_rate_role='AUTHORITATIVE_CASE_BOUNDARY_TARGET_NOT_MEASURED_REALIZED_FLUX',
 frozen_snapshot=f,chain=records)
(S/'validation/AUTHORITATIVE_INLET_FLOW_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
assert all(checks.values()),checks
local=json.loads((S/'provenance/INHERITED_INLET_FLUX_AUDIT.json').read_text())
(S/'validation/FROZEN_FIELD_LOCAL_CROSSCHECK.json').write_text(json.dumps(dict(role='FROZEN_FIELD_LOCAL_CROSSCHECK_ONLY',status='KNOWN_INPUT_LIMITATION',valid_area_percent=83.62354517646728,unsupported_area_percent=16.37645482353272,local_Q_positive_m3_s=local['quadrature'][-1]['partial_Q_positive_m3_s'],total_Q_source='AUTHORITATIVE_INLET_FLOW_AUDIT.json',no_fill_no_extrapolation=True),indent=2)+'\n')
oldgit=json.loads((S/'provenance/GIT_BEFORE.json').read_text());now={}
for repo in oldgit:
    def git(*a):return subprocess.check_output(['git','-C',repo,*a],env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'},text=True)
    now[repo]=dict(branch=git('branch','--show-current'),commit=git('rev-parse','HEAD'),status=git('status','--porcelain=v1','--untracked-files=all'))
assert now==oldgit,'STOP_BASELINE_DRIFT_GIT'
(S/'provenance/GIT_BEFORE.json').write_text(json.dumps(now,indent=2)+'\n')
for n in ['DONOR_HASHES_BEFORE.json','ORIGINAL_HASHES_BEFORE.json']:
    data=json.loads((S/'provenance'/n).read_text())
    for p,h in data.items():assert sha(Path(p))==(h['sha256'] if isinstance(h,dict) else h),(p,h)
print(json.dumps(dict(status='PASS',authoritative_Q=Q,checks=checks),indent=2))
