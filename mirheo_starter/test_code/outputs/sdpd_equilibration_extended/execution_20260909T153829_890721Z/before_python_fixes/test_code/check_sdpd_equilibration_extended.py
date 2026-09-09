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


def protection(audit):
    audit=Path(audit);before=read_json(audit/'source_before.json');modified=[];errors=[]
    for name,digest in before.items():
        p=PROJECT_ROOT/name
        if not p.is_file() or sha256_file(p)!=digest:
            (modified if name in ALLOWED_EXISTING_EDITS else errors).append(name)
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
    if errors:raise ValueError(str(result))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',default='test_code/outputs/sdpd_equilibration_extended/validation')
    p.add_argument('--protected-baseline',required=True)
    args=p.parse_args();destination=Path(args.output).resolve()
    if destination.exists():raise FileExistsError('CPU validation output must be a new immutable directory')
    destination.mkdir(parents=True)
    suite=unittest.defaultTestLoader.discover(str(PROJECT_ROOT/'test_code'),pattern='test_*.py',top_level_dir=str(PROJECT_ROOT))
    with (destination/'CPU_tests.log').open('x') as stream:
        result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if not result.wasSuccessful():raise RuntimeError('CPU_REGRESSION_FAILED_SEE_LOG')
    protected=protection(args.protected_baseline);write_json(destination/'protected_files.json',protected)
    from py_scripts.sdpd_diagnostics.equilibration import code_identity
    record={'status':'PASS','tests':result.testsRun,'code_sha256':code_identity(),
            'extended_config_sha256':sha256_file(PROJECT_ROOT/'py_scripts/sdpd_equilibration_extended.yaml'),
            'test_files_sha256':{str(f):sha256_file(f) for f in sorted((PROJECT_ROOT/'test_code').glob('test_*.py'))},
            'test_log_sha256':sha256_file(destination/'CPU_tests.log'),'protection':protected,
            'GPU_experiments':0,'real_restart_validation':'NOT_RUN','synthetic_CPU_data_is_experiment':False}
    write_json(destination/'CPU_validation.json',record)
    print(f'CPU PASS: {result.testsRun} tests; protected {protected["existing_files_checked"]} pre-existing files; {destination}')


if __name__=='__main__':main()
