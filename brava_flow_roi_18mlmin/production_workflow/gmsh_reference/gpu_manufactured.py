"""Execute frozen V2 through the existing formal BraVa GPU runner, without tuning."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
OLD = ROOT / 'manufactured_pipe/R1_GMSH_V2/cfd'
CASE = OLD.parent / 'cfd_gpu'
IDENTITY = ROOT / 'GPU_MANUFACTURED_RUN_IDENTITY.json'
STACK = Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def put(path, value):
    pending = Path(path).with_suffix(Path(path).suffix + '.pending')
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    pending.replace(path)


def science_options(text):
    tokens = shlex.split(text)
    for flag in ('-mat_type', '-vec_type'):
        index = tokens.index(flag)
        del tokens[index:index + 2]
    return [value for value in tokens if value != '-log_view_gpu_time']


def backend_evidence(text):
    # Only actual PETSc object views qualify; command-line requests do not.
    types = set(re.findall(r'^\s*type:\s*(\w+)\s*$', text, re.M))
    matrix = sorted(types & {'aijcusparse', 'seqaijcusparse', 'mpiaijcusparse'})
    vector = sorted(types & {'cuda', 'seqcuda', 'mpicuda'})
    return dict(verified=bool(matrix and vector), matrix_types=matrix,
                vector_types=vector, device_backend='CUDA' if matrix and vector else None)


def fatal_reason(text):
    if re.search(r'DIVERGED_\w+', text):
        return 'GPU_LINEAR_DIVERGENCE'
    if re.search(r'\b(?:nan|[+-]?inf)\b', text, re.I):
        return 'NONFINITE_NATIVE_OUTPUT'
    if re.search(r'CUDA error|CUDA_ERROR|out of memory|CUDART error|PETSC ERROR|'
                 r'error while loading shared libraries|lib.*not found', text, re.I):
        return 'GPU_BACKEND_OR_DEVICE_ERROR'
    return None


def prepare():
    assert not CASE.exists() and not IDENTITY.exists(), 'Never overwrite a run'
    authority = PROJECT / 'cases/balanced_pressure_final'
    production = read(authority / 'reports/execution.json')
    assert production['status'] == 'PASS'
    options = (authority / 'run/PETSC_OPTIONS.txt').read_text().strip()
    assert options == production['PETSC_OPTIONS']
    assert science_options(options) == science_options((OLD / 'run/PETSC_OPTIONS.txt').read_text())
    build = read(STACK / 'sv1_3q/reports/svmp_reuse_build.json')
    binary = PROJECT / 'gpu_solver_fix/svmultiphysics'
    library = Path(build['PETSc_build']['prefix']) / 'lib/libpetsc.so'
    assert sha(binary) == production['solver_sha256']
    assert sha(library) == production['PETSc_library_sha256']
    reuse = STACK / 'sv1_3q/external/svMultiPhysics-reuse/Code/Source/solver/sv13q_reuse.h'
    reuse_policy = STACK / 'sv1_3q/configs/policy.json'
    assert sha(reuse) == 'db635f95a7e0c187fa9768f1a836b57034bef4aca4319c2c24b0cd8fbddae220'
    assert sha(reuse_policy) == '410fb4d14dc817e8c96986d627ca322e8cc097557baf6109c88be000599da77a'
    hashes = read(OLD / 'input_hashes.json')
    for name, expected in hashes.items():
        assert sha(OLD / name) == expected, name
    qa = read(ROOT / 'GMSH_REFERENCE_MESH_V2_RESULT.json')
    assert qa['status'] == 'GMSH_REFERENCE_MESH_V2_PASS' and all(qa['checks'].values())
    assert hashes['SV_MESH/mesh-complete.mesh.vtu'] == qa['volume_sha256']
    assert sha(ROOT / 'gmsh_reference/candidate_v2/mesh/volume.vtu') == qa['volume_sha256']
    sys.path.insert(0, str(PROJECT / 'vendor'))
    from flow_solver_support.transient_checks import checkpoint_audit
    checkpoint = checkpoint_audit(OLD / 'run/8-procs/stFile_010.bin', 10, .01, ranks=8)
    CASE.mkdir()
    (CASE / 'run').mkdir()
    (CASE / 'reports').mkdir()
    shutil.copytree(OLD / 'SV_MESH', CASE / 'SV_MESH')
    for name in ('solver.xml', 'inlet_exact_velocity.dat', 'outlet_exact_traction.vtp'):
        shutil.copy2(OLD / 'run' / name, CASE / 'run' / name)
    (CASE / 'run/PETSC_OPTIONS.txt').write_text(options + '\n')
    policy = read(OLD / 'policy.json')
    assert policy['MPI_ranks'] == production['MPI_ranks'] == 8
    policy.update(maximum_wall_time_s=14400,
                  resource_plan='Formal BraVa GPU production: 8 MPI ranks share visible GPU 0; OMP=1; unchanged RA/ASM/local CPU LU')
    put(CASE / 'policy.json', policy)
    shutil.copy2(OLD / 'reports/FROZEN_SOLVER_IDENTITY.json', CASE / 'reports/SCIENTIFIC_BASELINE_IDENTITY.json')
    input_hashes = {name: sha(CASE / name) for name in hashes}
    put(CASE / 'input_hashes.json', input_hashes)
    assert {k: v for k, v in input_hashes.items() if k not in ('policy.json', 'run/PETSC_OPTIONS.txt')} == {
        k: v for k, v in hashes.items() if k not in ('policy.json', 'run/PETSC_OPTIONS.txt')}
    command = production['command']
    assert command == [build['PETSc_build']['candidate_wrapper'], build['PETSc_build']['launcher'],
                       '-n', '8', build['PETSc_build']['candidate_wrapper'], str(binary), 'solver.xml']
    gpu = subprocess.check_output(['nvidia-smi', '--query-gpu=name,driver_version,memory.total,uuid',
                                   '--format=csv,noheader'], text=True).strip()
    runtime = subprocess.check_output(['/usr/local/cuda/bin/nvcc', '--version'], text=True).strip()
    identity = dict(status='PREPARED_GPU_ONLY', VALIDATION_ONLY=True,
        execution_authority='cases/balanced_pressure_final/reports/execution.json',
        execution_authority_sha256=sha(authority / 'reports/execution.json'),
        scientific_authority='manufactured_pipe/R1_GMSH_V2/cfd',
        case=str(CASE.relative_to(ROOT)), mesh_sha256=qa['volume_sha256'],
        scientific_input_hashes={k: v for k, v in input_hashes.items() if k not in ('policy.json', 'run/PETSC_OPTIONS.txt')},
        solver_xml_byte_identical=True, PETSc_scientific_options_equal=True,
        solver_binary=str(binary), solver_sha256=sha(binary), PETSc_library=str(library),
        PETSc_library_sha256=sha(library), PETSC_OPTIONS=options,
        PETSC_OPTIONS_sha256=sha(CASE / 'run/PETSC_OPTIONS.txt'),
        command=command, MPI_ranks=8, OMP_NUM_THREADS=1, CUDA_VISIBLE_DEVICES='0',
        rank_GPU_mapping='All 8 historical production ranks share the single visible GPU 0',
        use_gpu_aware_mpi=False, host=gpu, CUDA_toolkit=runtime,
        runtime_library_path=build['runtime_library_path'],
        environment_wrapper_sha256=sha(command[0]), MPI_launcher_sha256=sha(command[1]),
        preconditioner_reuse=dict(policy='RA', option='-sv_pc_adaptive_rebuild true',
            source=str(reuse), source_sha256=sha(reuse), reference_policy=str(reuse_policy),
            reference_policy_sha256=sha(reuse_policy), max_age_steps=5, growth_factor=1.5,
            execution_settings_authority='Formal BraVa execution, not the older one-rank policy experiment'),
        initial_state='GPU_FRESH_STEP0', restart_assessment=dict(
            native_checkpoint_valid=True, checkpoint_sha256=checkpoint['sha256'],
            complete_integration_history=True, same_binary_mesh_physics_ranks=True,
            selected=False, reason='No saved rank-local node-order mapping is available for an independent restart correspondence check. User permits fresh GPU when compatibility remains uncertain; no conversion.'),
        dt_s=policy['dt_s'], nonlinear_tolerance=policy['nonlinear_tolerance'],
        steady_policy={k: policy[k] for k in ('save_interval_steps', 'steady_last_intervals',
            'velocity_change_limit', 'flow_change_limit', 'mass_limit')},
        maximum_wall_time_s=14400, native_checkpoint_cadence=10, saved_VTU_cadence=10,
        runtime_backend=dict(verified=False), automatic_CPU_fallback=False, automatic_BraVa=False)
    put(IDENTITY, identity)
    put(CASE / 'reports/FROZEN_SOLVER_IDENTITY.json', identity)
    print(json.dumps(dict(status=identity['status'],case=identity['case'],MPI_ranks=8,
                         initial_state=identity['initial_state']), indent=2), flush=True)


def run():
    import psutil
    sys.path.insert(0, str(ROOT / 'scripts'))
    from managed_job import stop_tree
    identity = read(IDENTITY)
    assert identity['status'] == 'PREPARED_GPU_ONLY'
    assert not (CASE / 'run/solver.log').exists()
    identity.update(status='GPU_MANUFACTURED_RUNNING', started_unix=time.time())
    put(IDENTITY, identity)
    command = [sys.executable, str(ROOT / 'scripts/run_pipe_cfd.py'), '--case', str(CASE),
               '--reference-root', str(PROJECT)]
    proc = subprocess.Popen(command, start_new_session=True)
    offset = 0
    observed = {}
    log_text = ''
    samples = []
    failure = None
    start = time.monotonic()
    prior_cpu = 0.
    prior_time = start
    try:
        while proc.poll() is None:
            now = time.monotonic()
            cpu = 0.
            rss = 0
            try:
                parent = psutil.Process(proc.pid)
                for child in [parent] + parent.children(recursive=True):
                    try:
                        observed[(child.pid, child.create_time())] = child
                        times = child.cpu_times()
                        cpu += times.user + times.system
                        rss += child.memory_info().rss
                    except psutil.Error:
                        pass
            except psutil.Error:
                pass
            path = CASE / 'run/solver.log'
            chunk = ''
            if path.exists():
                with path.open(errors='replace') as stream:
                    stream.seek(offset)
                    chunk = stream.read()
                    offset = stream.tell()
                log_text += chunk
            evidence = backend_evidence(log_text)
            identity['runtime_backend'] = evidence
            failure = fatal_reason(chunk)
            if failure == 'GPU_BACKEND_OR_DEVICE_ERROR' and not evidence['verified']:
                failure = 'GPU_BACKEND_UNAVAILABLE'
            if not evidence['verified'] and len(re.findall(r'^\s*NS\s+\d+-', log_text, re.M)) >= 2:
                failure = 'GPU_BACKEND_UNAVAILABLE'
            if now - start >= 14400:
                failure = 'GPU_CFD_TIME_LIMIT'
            gpu = subprocess.check_output(['nvidia-smi', '--query-gpu=utilization.gpu,memory.used',
                                          '--format=csv,noheader,nounits'], text=True).strip()
            utilization, memory = map(float, gpu.split(','))
            samples.append(dict(elapsed_s=now-start, GPU_utilization_percent=utilization,
                GPU_memory_MiB=memory, CPU_cores=max(0.,cpu-prior_cpu)/max(now-prior_time,1e-6),
                tree_RSS_GiB=rss/2**30))
            prior_cpu, prior_time = cpu, now
            identity.update(wall_time_s=now-start, last_resource_sample=samples[-1],
                peak_VRAM_MiB_sampled=max(s['GPU_memory_MiB'] for s in samples),
                preconditioner_reuse_observed=bool(re.search(r'SV13Q_BEGIN.*policy=RA.*reuse=1', log_text)))
            put(IDENTITY, identity)
            if failure:
                stop_tree(proc, observed)
                break
            time.sleep(5)
    finally:
        if proc.poll() is None:
            stop_tree(proc, observed)
    code = proc.wait()
    put(CASE / 'reports/GPU_RESOURCE_SAMPLES.json', samples)
    identity.update(wrapper_exit_code=code, wall_time_s=time.monotonic()-start,
                    failure_reason=failure, completed_unix=time.time())
    execution_path = CASE / 'reports/execution.json'
    if execution_path.exists():
        result = read(execution_path)
        identity['execution_sha256'] = sha(execution_path)
        identity['final_step'] = result.get('final_step')
        passed = result['status'] == 'PASS' and code == 0 and failure is None and identity['runtime_backend']['verified']
    else:
        passed = False
    identity['status'] = 'GPU_MANUFACTURED_STEADY_PASS' if passed else failure or 'GPU_MANUFACTURED_CFD_FAIL'
    put(IDENTITY, identity)
    print(json.dumps({k: identity[k] for k in ('status', 'wall_time_s', 'runtime_backend', 'failure_reason')},indent=2),flush=True)
    return 0 if passed else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'run'])
    stage = parser.parse_args().stage
    sys.exit(prepare() if stage == 'prepare' else run())
