"""Run the original finite-size cohort only after both flow movies are complete."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time,traceback

ROOT=Path(__file__).resolve().parents[1]
MB=ROOT/'microbubble'
PY=sys.executable
CUDA_PY='/venv/main/bin/python'

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def dump(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n')

def run(name,program,*args,python=PY,extra=None):
    state=dict(stage=name,status='RUNNING',started_unix=time.time())
    dump(MB/'data/SEQUENCE_PROGRESS.json',state)
    env=os.environ.copy()
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VTK_DEFAULT_OPENGL_WINDOW='vtkEGLRenderWindow')
    env.update(extra or {})
    command=[python,'-B',str(program),*args]
    with (MB/'logs'/f'{name}.log').open('a') as out:
        out.write(json.dumps({'command':command,'started':state['started_unix']})+'\n');out.flush()
        subprocess.run(command,env=env,stdout=out,stderr=subprocess.STDOUT,check=True,cwd=MB)
    state.update(status='COMPLETE',ended_unix=time.time())
    dump(MB/'data/SEQUENCE_PROGRESS.json',state)
    print(name,'COMPLETE',flush=True)

def main():
    if (MB/'data/USER_STOP.json').exists():
        print('CANCELLED_BY_USER: explicit user resumption is required before trajectory work.',flush=True)
        return
    (MB/'data').mkdir(parents=True,exist_ok=True)
    (MB/'logs').mkdir(parents=True,exist_ok=True)
    active=json.loads((ROOT/'reports/ACTIVE_FLOW.json').read_text())
    assert active['status']=='PASS'
    flow=ROOT/active['case']/'frozen_flow/steady_flow.vtu'
    assert sha(flow)==active['flow_sha256']
    for pose in ['0','15']:
        base=ROOT/'visualization'/('candidate_'+pose)
        media=json.loads((base/'MEDIA_VALIDATION.json').read_text())
        assert media['source_flow_sha256']==active['flow_sha256']
        assert len(media['media'])==5
        for record in media['media']:
            assert record['frames']==432 and sha(base/record['file'])==record['sha256']
    config=dict(source_root='/workspace/particle9a5_formal_trajectories_20260925T080115Z_dt1ms',
                flow_path=str(flow),flow_sha256=active['flow_sha256'],
                local_source_root='/home/lzy/projects/ulm_particle_formal_p9a5',
                physical_units='SI',dt_s=.0005,count=1500)
    if (MB/'config.json').exists():assert json.loads((MB/'config.json').read_text())==config
    else:dump(MB/'config.json',config)
    scripts=MB/'scripts'
    if not (MB/'data/cohort.json').exists():run('01_source',scripts/'campaign.py','prepare')
    if not (MB/'data/gpu_mesh_validation.json').exists():run('02_cuda_mesh',scripts/'gpu_validate.py','mesh',python=CUDA_PY)
    if not (MB/'data/pilot.json').exists():run('03_reference_pilot',scripts/'campaign.py','pilot')
    if not (MB/'data/native_kernel_parity.json').exists() or not json.loads((MB/'data/native_kernel_parity.json').read_text()).get('PASS'):
        run('04_native_parity',scripts/'test_native.py')
    run('05_exact_parity',scripts/'verify_native_bytes.py')
    run('06_cohort',scripts/'production_entry.py')
    run('07_cuda_tracks',scripts/'gpu_validate.py','tracks',python=CUDA_PY)
    run('08_summary',scripts/'summarize.py')
    # Preview and animation stages are launched after visual inspection by the operator.
    dump(MB/'data/TRAJECTORY_SEQUENCE_COMPLETE.json',dict(status='COMPLETE',next='Render candidate 0 and 15 after preview inspection'))

if __name__=='__main__':
    try:main()
    except BaseException:
        dump(MB/'data/SEQUENCE_FAILURE.json',dict(status='FAILED',traceback=traceback.format_exc(),time_unix=time.time()))
        raise
