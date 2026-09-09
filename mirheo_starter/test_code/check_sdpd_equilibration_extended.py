"""CPU-only validation certificate and protected-file audit; never launch CUDA."""
import argparse
import io
import os
from pathlib import Path
import subprocess
import unittest
from py_scripts.fluid_physics.common import PROJECT_ROOT,read_json,write_json,sha256_file,environment_identity

ALLOWED_EXISTING_EDITS={
    'py_scripts/sdpd_diagnostics/equilibration.py','py_scripts/sdpd_diagnostics/gpu_worker.py',
    'py_scripts/fluid_physics/budget_authorizations.py',
    'test_code/review_sdpd_equilibration.py','test_code/check_sdpd_diagnostics_browser.cjs'}


def execution_protection(audit):
    """Audit an explicitly approved pool append, while preserving all old data."""
    from py_scripts.fluid_physics.budget_authorizations import authorization_records
    audit=Path(audit);before=read_json(audit/'protected_before.json');errors=[]
    pool_path=PROJECT_ROOT/'runs/fluid_calibration/shared_budget_pool.json'
    old=read_json(audit/'shared_budget_pool.json.before');pool=read_json(pool_path)
    for key in set(old)|set(pool):
        if key not in ['members','authorization_records'] and old.get(key)!=pool.get(key):errors.append('POOL_FIELD_CHANGED '+key)
    for key in ['members','authorization_records']:
        if pool[key][:-1]!=old[key] or len(pool[key])!=len(old[key])+1:errors.append('POOL_NOT_SINGLE_APPEND '+key)
    records=authorization_records(pool);approval=read_json(audit/'approval_context.json')
    record=records[-1];request=read_json(Path(approval['approved_budget_request_file']))
    if (record['budget_request_sha256']!=approval['approved_budget_request_sha256'] or
            record['scope']!=request['authorization_scope'] or
            record['additional_seconds']!=request['additional_request_whole_seconds'] or
            pool['members'][-1]['directory']!=record['scope']['campaign_directory']):
        errors.append('APPENDED_AUTHORIZATION_NOT_APPROVED_SCOPE')
    for name,digest in before.items():
        path=Path(name)
        if path==pool_path:continue  # justified by the explicit append audit above
        if not path.is_file() or sha256_file(path)!=digest:errors.append('OLD_FILE_CHANGED '+name)
    completed=read_json(audit/'completed_GPU_raw_before_fixes.json')
    for name,digest in completed.items():
        if not (PROJECT_ROOT/name).is_file() or sha256_file(PROJECT_ROOT/name)!=digest:errors.append('GPU_EVIDENCE_CHANGED '+name)
    source=read_json(audit/'before_python_fixes/source_sha256.json');modified=[]
    allowed={'py_scripts/sdpd_diagnostics/extended.py','py_scripts/sdpd_diagnostics/extended_restart.py',
        'py_scripts/sdpd_diagnostics/extended_worker.py','test_code/review_sdpd_equilibration_extended.py',
        'test_code/check_sdpd_diagnostics_browser.cjs','test_code/check_sdpd_equilibration_extended.py'}
    for name,digest in source.items():
        if not (PROJECT_ROOT/name).is_file() or sha256_file(PROJECT_ROOT/name)!=digest:
            (modified if name in allowed else errors).append(name)
    if environment_identity()!=read_json(audit/'environment_before.json'):errors.append('NATIVE_ENVIRONMENT_CHANGED')
    result={'status':'FAIL' if errors else 'PASS','old_files_checked':len(before),
        'completed_GPU_files_checked':len(completed),'authorized_pool_append':not errors,
        'approved_additional_seconds':record['additional_seconds'],'modified_Python_viewer_files':modified,
        'old_raw_reports_ledgers_unchanged':not errors,'completed_GPU_raw_frozen_workers_unchanged':not errors,
        'unexpected_changes':errors,'analysis_launched_GPU':False}
    if errors:raise ValueError(str(result))
    return result


def protection(audit,execution_audit=None):
    audit=Path(audit);before=read_json(audit/'source_before.json');modified=[];errors=[]
    executed=execution_protection(execution_audit) if execution_audit else None
    for name,digest in before.items():
        p=PROJECT_ROOT/name
        if not p.is_file() or sha256_file(p)!=digest:
            allowed=name in ALLOWED_EXISTING_EDITS or (executed and name=='runs/fluid_calibration/shared_budget_pool.json')
            (modified if allowed else errors).append(name)
    expected=read_json(audit/'environment_before.json')
    if environment_identity()!=expected:errors.append('NATIVE_ENVIRONMENT_IDENTITY')
    git={}
    for name,path in [('native',PROJECT_ROOT/'vendor/Mirheo'),('vascular',PROJECT_ROOT.parent/'ulm_3D_vascular'),
                      ('sync',PROJECT_ROOT.parent/'ulm-3d-vascular-model-generation-sync-20260908T170007Z')]:
        git[name]=subprocess.check_output(['git','-c','filter.lfs.process=','-c','filter.lfs.required=false','-c','filter.lfs.clean=cat','-c','filter.lfs.smudge=cat',
            'status','--porcelain=v1','--untracked-files=all'],cwd=path,text=True,env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'})
    if git!=read_json(audit/'git_before.json'):errors.append('GIT_STATUS_CHANGED')
    result={'status':'FAIL' if errors else 'PASS','existing_files_checked':len(before),
            'authorized_minimal_existing_Python_viewer_changes':modified,'unexpected_changes':errors,
            'old_results_raw_reports_and_ledgers_unchanged':not errors,'native_environment_unchanged':environment_identity()==expected,
            'Git_status_unchanged':git==read_json(audit/'git_before.json'),'GPU_started':False}
    if executed:result['approved_execution_and_reanalysis_protection']=executed
    if errors:raise ValueError(str(result))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',default='test_code/outputs/sdpd_equilibration_extended/validation')
    p.add_argument('--protected-baseline',required=True)
    p.add_argument('--execution-baseline',help='Post-run audit directory with explicit approval and immutable before hashes; no GPU execution.')
    args=p.parse_args();destination=Path(args.output).resolve()
    if destination.exists():raise FileExistsError('CPU validation output must be a new immutable directory')
    destination.mkdir(parents=True)
    suite=unittest.defaultTestLoader.discover(str(PROJECT_ROOT/'test_code'),pattern='test_*.py',top_level_dir=str(PROJECT_ROOT))
    with (destination/'CPU_tests.log').open('x') as stream:
        result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if not result.wasSuccessful():raise RuntimeError('CPU_REGRESSION_FAILED_SEE_LOG')
    protected=protection(args.protected_baseline,args.execution_baseline);write_json(destination/'protected_files.json',protected)
    from py_scripts.sdpd_diagnostics.equilibration import code_identity
    record={'status':'PASS','tests':result.testsRun,'code_sha256':code_identity(),
            'extended_config_sha256':sha256_file(PROJECT_ROOT/'py_scripts/sdpd_equilibration_extended.yaml'),
            'test_files_sha256':{str(f):sha256_file(f) for f in sorted((PROJECT_ROOT/'test_code').glob('test_*.py'))},
            'test_log_sha256':sha256_file(destination/'CPU_tests.log'),'protection':protected,
            'GPU_experiments':0,'real_restart_validation':'EXISTING_REAL_FAILURE_NOT_RETRIED' if args.execution_baseline else 'NOT_RUN',
            'synthetic_CPU_data_is_experiment':False}
    write_json(destination/'CPU_validation.json',record)
    print(f'CPU PASS: {result.testsRun} tests; protected {protected["existing_files_checked"]} pre-existing files; {destination}')


if __name__=='__main__':main()
