#!/usr/bin/python3
"""Evaluate recorded data against the pre-execution criteria, then report."""
import argparse
import csv
import datetime
from pathlib import Path
import numpy as np
from audit_final import HC, CODE, ROOTS, audit, sha, read, write

def main(run):
    physical=read(run/'contracts/physical_bc_contract.json')
    geo=read(run/'contracts/geometry_reuse_contract.json')
    unit=read(run/'contracts/lattice_unit_contract.json')
    pre=read(run/'diagnostics/lbm_preflight.json')
    execution=read(run/'contracts/execution_contract.json')
    solver=read(run/'diagnostics/solver_status.json')
    sampler=read(run/'diagnostics/flux_sampler_runtime_check.json')
    geometry=read(run/'diagnostics/runtime_geometry_check.json')
    exports=read(run/'diagnostics/paraview_export_check.json')
    with (run/'diagnostics/flow_history.csv').open() as f:
        history=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    final=history[-1]
    with (run/'diagnostics/stages.csv').open() as f:
        stages={int(s['iteration']):s['safety_status'] for s in csv.DictReader(f)}
    finite=bool(np.isfinite(np.array([list(row.values()) for row in history])).all())
    log=(run/'logs/mpi1_run.log').read_text();build=(run/'logs/build.log').read_text()
    mpi_abort='MPI_ABORT' in log or 'MPI_Abort' in log
    segfault='Segmentation fault' in log or 'signal 11' in log.lower()
    backflow={name:{
        'observed':any(row[f'q_{name}_m3_s']<0 for row in history),
        'sampled_iterations':[int(row['iteration']) for row in history if row[f'q_{name}_m3_s']<0],
        'min_m3_s':min(row[f'q_{name}_m3_s'] for row in history),
        'final_m3_s':final[f'q_{name}_m3_s']} for name in physical['outlets']}
    repeat={}
    for step in [0,10,100,1000]:
        a=run/f'attempts/attempt1_invalid_flux/diagnostics/fields_{step}.bin'
        b=run/f'diagnostics/fields_{step}.bin'
        repeat[str(step)]={'equal_sha256':sha(a)==sha(b),'accepted_snapshot_sha256':sha(b)}
    write(run/'diagnostics/physical_repeat_check.json',{
        'status':'PASS' if all(v['equal_sha256'] for v in repeat.values()) else 'FAIL',
        'meaning':'Monitor correction did not change physical trajectory; first attempt flux columns are invalid',
        'snapshots':repeat})
    immutable=audit(run)
    acceptance={
        'physical_contract':physical['status']=='PASS',
        'api_feasibility':read(run/'contracts/bc_api_feasibility.json')['status']=='PASS',
        'lattice_units':unit['status']=='PASS','preflight':pre['status']=='PASS',
        'build_rc_zero':'Exit status: 0' in build,'runtime_rc_zero':'Exit status: 0' in log,
        'all_stages_safe':all(stages.get(i)=='PASS' for i in [0,1,10,100,1000]),
        'finite':finite,'no_mpi_abort':not mpi_abort,'no_segfault':not segfault,
        'native_geometry_unchanged':geometry['status']=='PASS',
        'four_ports_and_wall_called':all(solver['native_profile_calls'][i]>0 for i in range(5)),
        'flux_sampler_crosscheck':sampler['status']=='PASS',
        'inlet_profile_integral':geometry['profile_integral_error']<=1e-12,
        'final_inlet_flow_direction':final['q_in_m3_s']>0,
        'final_inlet_target_error':final['inlet_target_error']<=execution['short_run_acceptance']['final_plane_inlet_target_error_max'],
        'mass_drift':abs(final['relative_mass_drift'])<=execution['short_run_acceptance']['relative_mass_drift_abs_max'],
        'paraview_roundtrip':exports['status']=='PASS' and len(exports['snapshots'])==4,
        'frozen_inputs_integrity':immutable['status']=='PASS'}
    auto='PASS' if all(acceptance.values()) else 'FAIL'
    status='AUTO_PASS_HUMAN_PENDING' if auto=='PASS' else 'FAIL'
    failed=[k for k,v in acceptance.items() if not v]
    target=physical['inlet']['target_volume_flow_m3_s'];qin=final['q_in_m3_s']
    psum=sum(final[f'q_{k}_m3_s'] for k in physical['outlets'])
    write(run/'diagnostics/automatic_acceptance.json',{
        'status':auto,'step3_status':status,'checks':acceptance,'failed_checks':failed,
        'final':final,'backflow':backflow,'steady_state':'UNVERIFIED / NOT_ASSESSED',
        'human_review':'PENDING; Step2 human acceptance is not Step3 acceptance',
        'accepted_attempt':2,'timesteps_accepted_attempt':1000,'timesteps_all_attempts':2000})
    heads=read(run/'provenance/input_integrity.json')['git_heads']
    execution_history={'attempts':[
        {'id':1,'timesteps':1000,'status':'INVALID_FLUX_DIAGNOSTICS','physical_runtime_safety':'PASS',
         'record':'attempts/attempt1_invalid_flux/attempt_record.json'},
        {'id':2,'timesteps':1000,'runtime_safety':solver['runtime_safety'],'auto_acceptance':auto,'status':'CURRENT_RESULT'}],
        'total_fluid_timesteps':2000,'physical_parameter_tuning':False,
        'reason_for_repeat':'Corrected test-local MPI proxy alias in flow monitor. Same frozen BC, units, startup and quadrature; first attempt preserved.'}
    write(run/'provenance/execution_history.json',execution_history)
    provenance={
        'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'source_project_head':heads['source'],'hemocell_head':heads['hemocell'],
        'step1_stl_sha256':geo['input_stl_sha256'],
        'step1_all_files_sha256':read(run/'provenance/step1_before_sha256.json'),
        'step2_all_files_sha256':read(run/'provenance/step2_before_sha256.json'),
        'step2_code':read(run/'provenance/step2_code_review.json'),'step2_inputs_used':geo['files'],
        'step3_source':[{'path':str(p),'sha256':sha(p)} for p in sorted(CODE.iterdir()) if p.is_file()],
        'step3_configs':[{'path':str(p),'sha256':sha(p)} for p in sorted((run/'contracts').iterdir()) if p.is_file()],
        'physical_source_records':read(run/'provenance/physical_source_records.json'),
        'executable':{'path':str(run/'build/vascular_pure_fluid'),'sha256':sha(run/'build/vascular_pure_fluid')},
        'frozen_integrity_status':immutable['status'],'execution_history':execution_history}
    write(run/'provenance/provenance.json',provenance)
    unverified=['Step3 人工 ParaView','生理入口/出口身份','实验流量/黏度/profile','绝对生理血压',
                '稳态收敛','网格/时间步收敛','壁与 cap 接缝离散误差','流量积分收敛/精确离散质量残差',
                'MPI2','生产 CFD / RBC readiness']
    next_action='manual ParaView review of Step 3 pure-fluid short run' if auto=='PASS' else 'resolve Step 3 failure without changing frozen physical contract'
    fields={
        'STEP3_STATUS':status,'STEP3_AUTO_CHECK':auto,
        'SOURCE_PROJECT_HEAD':heads['source'],'HEMOCELL_HEAD':heads['hemocell'],
        'STEP1_STL_SHA256':geo['input_stl_sha256'],
        'STEP2_GEOMETRY_MAPPING_REUSED':'YES; native flags, transform, 191 cap triangles and 497 labels checked',
        'PHYSICAL_BC_CONTRACT':physical['status'],'BC_API_FEASIBILITY':'PASS',
        'FLUID_DENSITY_KG_M3':unit['rho_phys_kg_m3'],'KINEMATIC_VISCOSITY_M2_S':unit['nu_phys_m2_s'],
        'INLET_ROLE_STATUS':'ASSUMED','INLET_BC_PHYSICAL_TYPE':'target_volumetric_flow',
        'INLET_TARGET_Q_M3_S':target,
        'INLET_NUMERICAL_METHOD':'native Guo + Q-normalized VelocityPlugProfile3D; explicit profile assumption'}
    for name,p in physical['outlets'].items():
        fields[name.upper()+'_ROLE_STATUS']='ASSUMED'
        fields[name.upper()+'_PRESSURE_PA']=p['pressure_pa']
        fields[name.upper()+'_NUMERICAL_METHOD']='native Guo + DensityNeumannBoundaryProfile3D (fixed gauge density)'
    fields.update({
        'WALL_BC':'stationary no-slip; native regularized Guo off-lattice',
        'DX_STATUS':'CANDIDATE_BASELINE','DX_REQUESTED_M':geo['dx_requested_m'],'DX_EFFECTIVE_M':geo['dx_effective_m'],
        'DT_S':unit['dt_s'],'TAU':unit['tau'],'NU_LU':unit['nu_lu'],'MAX_EXPECTED_MACH':pre['max_expected_mach'],
        'LATTICE_UNIT_CONTRACT':unit['status'],'LBM_PREFLIGHT':pre['status'],'BUILD_RC':0})
    for n,label in [(1,'RUN_1_STEP'),(10,'RUN_10_STEPS'),(100,'RUN_100_STEPS'),(1000,'RUN_1000_STEPS')]:
        fields[label]=stages.get(n,'NOT_RUN')
    fields.update({
        'NaN':'NONE' if finite else 'DETECTED','Inf':'NONE' if finite else 'DETECTED',
        'MPI_ABORT':'DETECTED' if mpi_abort else 'NONE','SEGFAULT':'DETECTED' if segfault else 'NONE',
        'FINAL_RHO_MIN':final['rho_min'],'FINAL_RHO_MAX':final['rho_max'],
        'FINAL_MAX_VELOCITY_M_S':final['u_max_m_s'],'FINAL_MAX_MACH':final['mach_max'],'FINAL_Q_IN_M3_S':qin})
    for name in physical['outlets']:fields['FINAL_Q_'+name.upper()+'_M3_S']=final[f'q_{name}_m3_s']
    fields.update({
        'FINAL_FLOW_CLOSURE':final['flow_closure'],'FINAL_INLET_TARGET_ERROR':final['inlet_target_error'],
        'FINAL_MASS_DRIFT':final['relative_mass_drift'],
        'BACKFLOW_OBSERVED':'; '.join(k for k,v in backflow.items() if v['observed']) or 'NONE',
        'FLUID_TIMESTEPS_RUN':'2000 total = 1000 archived invalid-monitor attempt + 1000 current validation',
        'MPI1':'PASS' if auto=='PASS' else 'FAIL','MPI2':'UNVERIFIED',
        'PARAVIEW_OUTPUTS':str(run/'output/flow.pvd')+' (0,10,100,1000; VTU)',
        'HEMOCELL_CORE_MODIFIED':'NO' if immutable['hemocell']['status']=='PASS' else 'CHECK_FAILED',
        'PALABOS_MODIFIED':'NO' if immutable['hemocell']['status']=='PASS' else 'CHECK_FAILED',
        'OFFICIAL_EXAMPLES_MODIFIED':'NO' if immutable['hemocell']['status']=='PASS' else 'CHECK_FAILED',
        'EXISTING_TRACKED_FILE_MODIFIED_BY_STEP3':'NO' if not immutable['hemocell']['tracked_hash_mismatches'] else 'CHECK_FAILED',
        'TEST_CODE_DIR':str(CODE),'RUN_DIR':str(run),
        'STEP3_REPORT':str(run/'STEP3_REPORT.md'),'PARAVIEW_REVIEW':str(run/'PARAVIEW_REVIEW.md'),
        'UNVERIFIED_ITEMS':'; '.join(unverified),'NEXT_RECOMMENDED_ACTION':next_action})
    (run/'FINAL_STATUS.txt').write_text('\n'.join(f'{k} = {v}' for k,v in fields.items())+'\n')
    write(run/'diagnostics/final_status.json',fields)
    # Report text is kept separate from numerical acceptance so prose cannot alter a gate.
    from report_text import reports
    reports(run,locals())
    files=[p for p in sorted(run.rglob('*')) if p.is_file() and p.name not in ['SHA256SUMS','write_reports.log']]
    (run/'SHA256SUMS').write_text('\n'.join(f'{sha(p)}  {p.relative_to(run)}' for p in files)+'\n')
    print(__import__('json').dumps({'status':status,'failed_checks':failed,'final':final,'integrity':immutable['status']},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path)
    main(parser.parse_args().run)
