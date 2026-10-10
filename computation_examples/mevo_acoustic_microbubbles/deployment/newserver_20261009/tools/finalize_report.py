"""Collect migration evidence; fail if any validation phase is incomplete."""
from pathlib import Path
import datetime
import hashlib
import json
import subprocess

r=Path('/workspace/mevo_acoustic_microbubbles')
m=r/'migration'
def read(path):
    return json.loads(path.read_text())
v=read(m/'validation/VALIDATION.json')
assert v['status']=='PASS', v
assert len(v['phases'])==7 and all(x['status']=='PASS' for x in v['phases'].values())
files=read(m/'inventory/SOURCE_FILES.json')
refs=read(m/'inventory/REFERENCE_FILES.json')
for row in files:
    assert hashlib.sha256((r/'repo'/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
for row in refs:
    assert hashlib.sha256((r/'references'/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
git=lambda *args: subprocess.check_output(['git','-C',str(r/'repo'),*args],text=True).strip()
assert git('rev-parse','HEAD')=='66d4cebec18555a26bc9a8035cc3082756d0e76f'
assert git('status','--porcelain')=='',git('status','--porcelain')
gpu=read(m/'validation/f5_gpu_checks.json')['detail']
flow=read(m/'validation/f5'/gpu['label']/'FLOW_GPU_VALIDATION.json')
saved=read(m/'validation/saved_data.json')['detail']
render=read(m/'validation/render_video.json')['detail']
report=dict(
    status='COMPLETE_VERIFIED', completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    endpoint='root@95.3.33.46:15283', local_ssh_alias='mevo-new', root=str(r),
    project=str(r/'project'), source_branch=git('branch','--show-current'),
    source_commit=git('rev-parse','HEAD'), source_worktree_clean=True,
    original_science_commit=v['original_science_commit'],
    data=dict(source_files=len(files),source_logical_bytes=sum(f['size'] for f in files),
              reference_files=len(refs),reference_bytes=sum(f['size'] for f in refs)),
    host=read(m/'inventory/HOST.json'),
    compute_environment=str(r/'.venv'),render_environment=str(r/'.venv-render'),
    validation=v, f5_check_groups=len(gpu['checks']), f5_check_label=gpu['label'],
    original_assertions_unchanged=gpu['original_assertions_unchanged'],
    original_ON_OFF_prefixes_bitwise_identical=gpu['archived_ON_OFF_prefixes_bitwise'],
    max_GPU_CPU_velocity_error_m_s=max(x['max_velocity_abs_error_m_s'] for x in flow['groups'].values()),
    hdf5_file_paths_fully_read=saved['hdf5_files'],hdf5_datasets_read=saved['hdf5_datasets'],
    saved_trajectory_CPU_points=saved['CPU_actual_saved_flow_points'],
    saved_trajectory_CPU_error_m_s=saved['CPU_flow_error_m_s'],
    offscreen_render_passed=render['offscreen_render'],
    original_1201_frame_video_decode_passed=render['full_decode'],
    h264_encode_passed=render['h264_encode'],
    resolutions=[
      'Default Warp 1.18.0 CUDA13.4 wheel could not use driver 550.120; official 1.18.0+cu12 wheel verified and tested',
      'Original renderer needs newer PyVista API; separate original-version render environment created; compute dependencies unchanged',
      'F1/F4 missing input references restored with relative SHA-identical links; absolute historical provenance left untouched'],
    science_results_modified=False,historical_jobs_resumed=False,new_full_scientific_runs=0,
    global_environment_modified=False,host_driver_modified=False,
    scope=read(m/'inventory/MIGRATION_SCOPE.json'),
    local_report='/home/lzy/projects/server_migration/mevo_newserver_20261009/MIGRATION_REPORT.json',
    local_complete_backup='/home/lzy/projects/server_migration/vast4090_20261009')
(r/'MIGRATION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['status','f5_check_groups','max_GPU_CPU_velocity_error_m_s','hdf5_file_paths_fully_read','offscreen_render_passed']}))
