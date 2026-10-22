"""Close the explicit geometry stop gate; summarize completed supervision and identities."""
from common import *
import datetime
import os
import psutil
from PIL import Image

def main():
    assert read(HERE/'TEST_RESULTS.json')['status']=='PASS'
    assert read(HERE/'P2_CURVED_GEOMETRY_AUDIT.json')['status']=='P2_CURVED_GEOMETRY_NOT_SUPPORTED'
    assert not (HERE/'TASK.json').exists()
    active=[]
    for process in psutil.process_iter(['pid','cmdline']):
        args=process.info['cmdline'] or []
        if process.pid!=os.getpid() and any(str(HERE) in a and a.endswith('/supervise.py') for a in args):active.append(process.pid)
    assert not active,active
    # Final packaging checks after any plot-layout-only changes.
    for f in read(HERE/'FIGURE_MANIFEST.json')['figures']:
        assert sha(HERE/f['file'])==f['sha256']
        with Image.open(HERE/f['file']) as im:assert abs(im.info['dpi'][0]-300)<.1
    for name,digest in read(HERE/'INPUT_SHA256.json').items():
        if not (BASE/name).is_relative_to(HERE):assert sha(BASE/name)==digest,name
    gate=dict(status='P2_CURVED_GEOMETRY_NOT_SUPPORTED',
        P2_P1_NATIVE_HIGHER_ORDER_WORTH_CONTINUING='NOT_EVALUABLE',
        P2_REPRESENTATION_ADVANTAGE='SUPPORTED',
        final_recommendation='P2P1_NOT_EVALUABLE_DUE_TO_IMPLEMENTATION_BLOCKER',
        stopped_at='D2 / Phase E 9.1 prelaunch geometry gate',
        phases={'A':'NATIVE_P2P1_PATH_CONFIRMED','B':'PASS','C':'PASS','D1':'PASS','D2':'BLOCKED',
                'E':'NOT_RUN','F':'NOT_RUN'},
        reason='Exact retained native curved assembly uses different pressure/velocity mappings and overwrites the first continuity quadrature Jacobian with the affine pressure Jacobian. Direct construct_fluid residual captures reproduce this inconsistency in two valid elevated V2 cells.',
        stop_is_not='Not a failed GPU run, not proof that P2 approximation lacks value, and not evidence that the original P1 anomaly is caused by this unexecuted P2 path.',
        source_patch=False,pipe_CFD_runs=0,BraVa_CFD_runs=0,particle_runs=0,full_curved_meshes=0,
        automatic_followup=False,production_certified=False,trajectory_fixed=False)
    put('P2_DECISION_GATE.json',gate)
    for name in ['P2_PIPE_FIELD_VALIDATION.json','P2_PIPE_RUN_IDENTITY.json','P2_BRAVA_FIELD_VALIDATION.json','P2_BRAVA_RUN_IDENTITY.json']:
        put(name,dict(status='NOT_RUN',reason=gate['reason'],blocked_at=gate['stopped_at'],binary_executed_as_CFD=False,
            new_field=None,new_checkpoint=None,wall_time_s=None,GPU_memory=None,MPI_layout=None,steady_evidence=None))
    records=[]
    for p in sorted((HERE/'monitoring').glob('*/EXECUTION.json')):
        r=read(p)
        records.append(dict(name=p.parent.name,receipt_path=str(p.relative_to(HERE)),receipt_sha256=sha(p),
            exit_code=r['exit_code'],stop_reason=r['stop_reason'],elapsed_s=r['elapsed_s'],
            sampled_peak_tree_RSS_GiB=r['peak_tree_RSS_GiB'],completed_UTC=r['completed_UTC'],
            samples=len(r['samples']),initial_memory=r['initial_memory']))
    receipt=dict(status='COMPLETE_AT_EXPLICIT_STOP_GATE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        jobs=records,job_count=len(records),cumulative_supervised_non_CFD_s=sum(r['elapsed_s'] for r in records),
        non_CFD_budget_s=5400,sampled_peak_tree_RSS_GiB=max(r['sampled_peak_tree_RSS_GiB'] for r in records),
        RSS_sampling_period_s=1,own_active_supervisors=active,TASK_json_present=False,automatic_next_job=False,
        pipe_CFD_runs=0,BraVa_CFD_runs=0,particle_runs=0,GPU_CFD_or_particle_work='NOT_RUN',
        CPU_layout='Single-process NumPy/geometry/plots; native unit programs initialize one MPI rank, no CFD or linear solve; BLAS/OMP threads1.',
        initial_harness_failures=[dict(job='native_build',reason='Test harness enum name corrected; production unchanged.'),
            dict(job='native_build_attempt2',reason='Test harness initially used reduced-pressure quadrature initializer for TET10. Corrected to actual native select_ele -> init_fs_msh chain; production unchanged.')],
        short_operations='Read-only preflight, source inspection, transfers and final JSON/Markdown packaging are not represented as scientific supervised runtime.',
        final_packaging_checks={'figure_SHA_and_300dpi':'PASS','protected_historical_SHA':'PASS'},
        others_services_modified=False,system_libraries_modified=False,cgroup_modified=False)
    put('EXECUTION_RECEIPT.json',receipt)
    source=read(HERE/'SOURCE_IDENTITY.json');pipe=read(HERE/'P2_POISEUILLE_REPRESENTATION.json')
    identity=dict(status='VALIDATION_ONLY_IMPLEMENTATION_BLOCKED',local_task_start_HEAD=read(HERE/'PROTOCOL.json')['production_source_commit_at_task_start'],
        native_source_commit=source['source_commit'],native_binary_sha256=source['binary_sha256'],PETSc=source['PETSc'],
        current_CUDA_library_identity='CUDA_LIBRARY_IDENTITY.json',
        native_harness='NATIVE_FIXTURE_BUILD.json',native_curved_harness='NATIVE_GEOMETRY_BUILD.json',
        mesh='Deterministic field-only affine elevation, same174785cells/36211corners; no full curved mesh',
        elevated_mesh_sha256=pipe['mesh']['elevation_sha256'],affine_runtime_artifact=pipe['mesh']['elevation_artifact'],
        input_fields='Analytic quadratic, near-wall polynomial and Poiseuille; existing M0/P1 fields read only; no new CFD output',
        original_geometry_BC_protocols_unchanged=True,protected_input_manifest='INPUT_SHA256.json',
        figures='FIGURE_MANIFEST.json',tests='TEST_RESULTS.json',gate='P2_DECISION_GATE.json',
        execution='EXECUTION_RECEIPT.json',validation_sampler_source_sha256=sha(HERE/'p2_sampler.py'),
        output_checkpoint=None,actual_CFD_MPI_layout=None,actual_CFD_steady_evidence=None,
        commit_identity='git log -- this directory; manifest does not include its own enclosing commit hash',push_performed=False)
    put('RESULT_IDENTITY.json',identity)
    # Input/source provenance includes all task programs; avoids circular output hashes.
    inputs=read(HERE/'INPUT_SHA256.json')
    for p in list(HERE.glob('*.py'))+list(HERE.glob('*.cpp'))+[HERE/'PROTOCOL.json']:
        inputs[str(p.relative_to(BASE))]=sha(p)
    put('INPUT_SHA256.json',inputs)
    print(gate['final_recommendation'],receipt['job_count'],receipt['cumulative_supervised_non_CFD_s'])

if __name__=='__main__':main()
