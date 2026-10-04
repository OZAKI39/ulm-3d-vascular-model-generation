"""Certify only a GPU-produced, natively saved state after frozen steady gates."""
from pathlib import Path
import json
import re
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT / 'scripts'), str(ROOT.parent / 'vendor')]
from gpu_manufactured import CASE, IDENTITY, read, put, sha
from pipe_measurements import PipeMeasurements
from flow_solver_support.steady_monitor import stop_gate
from flow_solver_support.transient_checks import checkpoint_audit
from flow_solver_support.solver_checks import linear_gate, nonlinear_gate
from validate_field import analyze


def main():
    identity = read(IDENTITY)
    assert identity['status'] == 'GPU_MANUFACTURED_STEADY_PASS'
    assert identity['runtime_backend']['verified']
    assert identity['preconditioner_reuse_observed']
    execution_path = CASE / 'reports/execution.json'
    execution = read(execution_path)
    assert sha(execution_path) == identity['execution_sha256']
    assert execution['solver_sha256'] == identity['solver_sha256']
    assert execution['PETSC_OPTIONS'] == identity['PETSC_OPTIONS']
    native_log = (CASE / 'run/solver.log').read_text(errors='replace')
    tolerances = re.search(r'tolerances: relative=([^,]+), absolute=([^,]+),', native_log)
    assert tolerances and tuple(map(float, tolerances.groups())) == (1e-10, 1e-24)
    identity['native_KSP_tolerances'] = dict(relative=1e-10, absolute=1e-24,
        source='First global KSP object view in native solver.log',
        log_sha256=sha(CASE / 'run/solver.log'))
    linear_gate(execution)
    policy = read(CASE / 'policy.json')
    nonlinear_gate(execution['history'], max_iterations=policy['maximum_nonlinear_iterations'],
                   tolerance=policy['nonlinear_tolerance'])
    assert stop_gate(execution['intervals'], execution['states'][-1], 0, 0, policy)
    final_path = CASE / execution['final_vtu']
    checkpoint = final_path.with_name('stFile_%03d.bin' % execution['final_step'])
    cp = checkpoint_audit(checkpoint, execution['final_step'], policy['dt_s'], ranks=8)
    measure = PipeMeasurements(CASE / 'SV_MESH/mesh_arrays.npz', policy['Q_target_m3_s'], policy['Umean_m_s'])
    u, p = measure.read(final_path)
    measured = measure.measure(u, p)
    assert stop_gate(execution['intervals'], measured, 0, 0, policy)
    steady = dict(final_step=execution['final_step'], final_VTU=str(final_path.relative_to(ROOT)),
        final_VTU_sha256=sha(final_path), final_native_checkpoint=str(checkpoint.relative_to(ROOT)),
        final_native_checkpoint_sha256=cp['sha256'], measured=measured,
        qualifying_saved_step=execution['stop']['step'],
        qualifying_five_intervals=execution['intervals'][-5:],
        velocity_pressure_fields_modified=False, CPU_step10_used_for_certification=False)
    identity['steady_saved_state'] = steady
    libraries_path = CASE / 'reports/GPU_LOADED_LIBRARIES.json'
    samples_path = CASE / 'reports/GPU_RESOURCE_SAMPLES.json'
    samples = read(samples_path)
    identity['loaded_GPU_libraries'] = read(libraries_path)
    identity['loaded_GPU_libraries_sha256'] = sha(libraries_path)
    identity['resource_record'] = dict(
        samples_sha256=sha(samples_path), sampling_period_s=5,
        GPU_memory_scope='Sampled total device memory, including existing non-job allocations; not exclusive process VRAM',
        GPU_utilization_scope='Sampled whole-device utilization; not a performance benchmark',
        CPU_scope='Sum of CPU time increments in this execution process tree, expressed as occupied cores',
        peak_device_VRAM_MiB=max(s['GPU_memory_MiB'] for s in samples),
        mean_GPU_utilization_percent=sum(s['GPU_utilization_percent'] for s in samples)/len(samples),
        mean_CPU_cores_sampled=sum(s['CPU_cores'] for s in samples)/len(samples),
        peak_tree_RSS_GiB=read(ROOT / 'monitoring/manufactured_gmsh_v2_gpu/EXECUTION.json')['peak_tree_RSS_GiB'])
    put(IDENTITY, identity)
    analyze(CASE, 'manufactured_gmsh_v2_gpu')
    result = read(ROOT / 'FIELD_VALIDATION.json')
    passed = result['status'] == 'REFERENCE_PIPE_FIELD_PASS'
    result.update(status='REFERENCE_PIPE_FIELD_PASS' if passed else 'REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC',
        GPU_manufactured_status='GPU_MANUFACTURED_STEADY_PASS',
        local_field_accuracy_assessed=True, field_certified=passed,
        GPU_execution_identity_sha256=sha(IDENTITY), steady_saved_state=steady,
        failure_stage=None if passed else 'SCIENTIFIC_FIELD_GATES',
        particle_cases_executed=0)
    for source in [Path(__file__), ROOT / 'manufactured_pipe/PROTOCOL.json',
                   HERE / 'candidate_v2/QA_SAMPLES.npz',
                   *[ROOT / 'scripts' / name for name in
                     ('pipe_measurements.py', 'pipe_probe_geometry.py', 'particle_field_inputs.py', 'poiseuille.py')]]:
        result['source_sha256'][str(source.relative_to(ROOT))] = sha(source)
    put(ROOT / 'FIELD_VALIDATION.json', result)
    print(json.dumps(dict(status=result['status'], checks=result['checks']), indent=2), flush=True)


if __name__ == '__main__':
    main()
