#!/usr/bin/env python3
"""Read-only source export into two independent, auditable publication snapshots."""
from pathlib import Path
import collections, gzip, hashlib, json, os, shutil, subprocess, time

BASE = Path(__file__).resolve().parent
SIM_REPO = Path('/home/lzy/projects/computation_examples/_archives/mevo_capture_f10ag_20261019')
SIM = SIM_REPO / 'computation_examples/mevo_acoustic_microbubbles'
VISION = Path('/home/lzy/projects/external_connect/vision_tracking_research')
MEDIA = {'.mp4','.avi','.mov','.mkv','.webm','.gif'}
ENV = {'.git','__pycache__','.pytest_cache','.mypy_cache','.ruff_cache','node_modules','_envs','.cache'}
MODEL_KEEP = {'20261016T202356_1416085f8f','20261019T212720_9f1dc0b82f','20261019T213744_3e56114286'}

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return h.hexdigest()

def write(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')

def git_identity(p):
    def run(*args):
        q=subprocess.run(['git','-C',str(p),*args],capture_output=True,text=True)
        return q.stdout.strip() if q.returncode==0 else None
    return dict(path=str(p),resolved=str(p.resolve()),git_root=run('rev-parse','--show-toplevel'),
                commit=run('rev-parse','HEAD'),branch=run('branch','--show-current'),
                status=run('status','--porcelain=v1','--untracked-files=all'))

def reason(kind,p,size):
    s=p.as_posix(); parts=p.parts
    if p.name in {'.DS_Store','Thumbs.db'} or p.suffix in {'.pyc','.pyo','.so','.dll','.exe','.whl','.o'}: return 'generated_or_platform_binary'
    if p.name in {'.env','id_rsa','id_ed25519','hosts.yml'} or p.suffix in {'.pem','.key'}: return 'private_local_configuration'
    if p.name in {'observation_read_request.json','observation_status.json'} and len(parts)==1: return 'transient_hardware_mailbox_archived_separately'
    if p.suffix in MEDIA: return 'full_video_or_animation_external_artifact; scientific_data_and_renderer_retained'
    if 'tools' in parts and kind=='chronos':return 'redistributable_VLC_binaries'
    if 'drivers' in parts and kind=='afg':return 'vendor_driver_installers'
    if kind=='sim':
        if '/lf/' in '/'+s and p.name=='trajectory.h5' and size>8_000_000:
            return 'bulky_low_fidelity_diagnostic_trajectory; HF_dataset_LF_summary_and_step_records_retained'
    if kind=='vision':
        if s.startswith('smart_online_tracking_v2/models/') or '/microscope_20261019_v1/models/' in s:
            return 'downloadable_base_model_or_recomputable_embedding; identities_retained'
        if s.startswith('third_party/checkpoints/') or s.startswith('collective_motion_v1/models/SEA-RAFT-S/') and size>10_000_000:
            return 'downloadable_upstream_model'
        if s.startswith('collective_motion_v2/data/') and p.suffix=='.npy' and size>1_000_000:
            return 'regenerable_synthetic_training_array; seeds_generator_labels_and_results_retained'
        if '/sam_auto_training/' in s and p.name=='decoder.pt' and not MODEL_KEEP.intersection(parts):
            return 'superseded_unpromoted_decoder; all_training_validation_and_promotion_logs_retained'
        if '/third_party/' in '/'+s and p.suffix in {'.png','.jpg','.jpeg','.bmp','.tif','.tiff','.npy','.npz','.ipynb'} and size>500_000:
            return 'upstream_demo_media_or_benchmark_asset'
        if s.startswith('collective_motion_v4/results/') and p.name=='records.jsonl.gz' and size>15_000_000:
            if not any(x in s for x in ['replay_final_v6/','full_video_verified_v5/']):
                return 'superseded_bulk_PIV_records; final_verified_records_and_all_summaries_retained'
    return None

class Export:
    def __init__(self,name):
        self.root=BASE/name;self.root.mkdir(exist_ok=True)
        self.files=[];self.excluded=[];self.aliases=[];self.identities=[];self.pruned=[]
    def tree(self,source,dest,kind):
        source=Path(source);self.identities.append(git_identity(source))
        for d,dirs,files in os.walk(source,followlinks=True):
            dp=Path(d);rel=dp.relative_to(source)
            for x in list(dirs):
                q=dp/x;rp=rel/x
                if x in ENV or x.startswith('.venv') or (kind=='chronos' and x=='tools') or (kind=='afg' and x=='drivers'):
                    dirs.remove(x);self.pruned.append(dict(source=str(q),reason='environment_cache_vcs_or_vendor_binary_tree'));continue
                if q.is_symlink():
                    if x=='models' and 'microscope_sessions' in rp.parts:
                        dirs.remove(x);target=os.path.relpath(self.root/dest/'smart_online_tracking_v2/models',self.root/dest/rel)
                        out=self.root/dest/rp;out.parent.mkdir(parents=True,exist_ok=True)
                        if not out.exists() and not out.is_symlink():out.symlink_to(target)
                        self.aliases.append(dict(path=str(Path(dest)/rp),source=str(q),source_target=os.readlink(q),portable_target=target));continue
                    self.aliases.append(dict(path=str(Path(dest)/rp),source=str(q),source_target=os.readlink(q),action='materialized_directory'))
            for name in sorted(files):
                p=dp/name;rp=rel/name
                if not p.is_file():continue
                self.file(p,Path(dest)/rp,kind,rp)
    def file(self,p,dest,kind,rel=None,force=False):
        p=Path(p);rel=rel or Path(dest);st=p.stat();why=None if force else reason(kind,rel,st.st_size)
        rec=dict(source=str(p),path=str(dest),source_bytes=st.st_size,source_mtime_ns=st.st_mtime_ns)
        if why:
            rec['reason']=why
            # Streaming hashing even for omitted videos; no frame decoding or RAM loading.
            rec['sha256']=sha(p);self.excluded.append(rec);return
        out=self.root/dest;out.parent.mkdir(parents=True,exist_ok=True)
        compress=(st.st_size>2_000_000 and p.suffix in {'.json','.jsonl','.csv','.log','.h5','.raw','.tar','.npy'})
        if compress:
            out=Path(str(out)+'.gz');rec['path']=str(out.relative_to(self.root));rec['encoding']='gzip_exact_original_bytes'
            h=hashlib.sha256()
            with p.open('rb') as i,out.open('wb') as o,gzip.GzipFile(filename='',mode='wb',fileobj=o,mtime=0,compresslevel=3) as g:
                for b in iter(lambda:i.read(4*1024*1024),b''):h.update(b);g.write(b)
            rec['source_sha256']=h.hexdigest()
        else:
            shutil.copy2(p,out);rec['encoding']='identity';rec['source_sha256']=sha(out)
        if p.stat().st_size!=st.st_size or p.stat().st_mtime_ns!=st.st_mtime_ns:raise RuntimeError('Source changed while copying: '+str(p))
        if out.stat().st_size>=95_000_000:raise RuntimeError('Review oversized artifact: '+str(out))
        rec['bytes']=out.stat().st_size;rec['sha256']=sha(out)
        rec['source_mode']=oct(st.st_mode & 0o777)
        self.files.append(rec)
    def finish(self):
        sync=self.root/'sync/20261022'
        write(sync/'SOURCE_IDENTITIES.json',self.identities)
        write(sync/'COPIED_FILES.json',self.files)
        write(sync/'EXCLUDED_FILES.json',self.excluded)
        write(sync/'PRUNED_DIRECTORIES.json',self.pruned)
        write(sync/'SYMLINK_POLICY.json',self.aliases)
        write(sync/'PAYLOAD_SUMMARY.json',dict(files=len(self.files),payload_bytes=sum(x['bytes'] for x in self.files),
            original_included_bytes=sum(x['source_bytes'] for x in self.files),excluded_files=len(self.excluded),
            excluded_bytes=sum(x['source_bytes'] for x in self.excluded),compressed_files=sum(x['encoding']!='identity' for x in self.files),
            largest_file_bytes=max(x['bytes'] for x in self.files)))
        shutil.copy2(__file__,sync/'export_snapshots.py')
        print(self.root.name, (sync/'PAYLOAD_SUMMARY.json').read_text(),flush=True)

def main():
    start=time.monotonic()
    sim=Export('simulation')
    sim.tree(SIM,'computation_examples/mevo_acoustic_microbubbles','sim')
    sim.tree(BASE/'audit/f10ag_server_logs','computation_examples/mevo_acoustic_microbubbles/sync/20261022/f10ag_server_logs','sim')
    sim.file('/home/lzy/projects/microbubble_aggregation.md','docs/shared_project_overview_20261022.md','shared',force=True)
    sim.finish()
    exp=Export('experiment')
    exp.tree(VISION,'external_connect/vision_tracking_research','vision')
    exp.tree('/mnt/f/CHRONOS_control','external_connect/CHRONOS_control','chronos')
    exp.tree('/mnt/f/Tektronix_control','external_connect/Tektronix_control','afg')
    # Keep the distinct WSL event log and the stopped-controller state as evidence only.
    for p in ['/home/lzy/projects/external_connect/Tektronix_control/control.log','/mnt/f/Tektronix_control/observation_status.json','/mnt/f/Tektronix_control/observation_read_request.json']:
        label='wsl_'+Path(p).name if p.startswith('/home') else 'windows_'+Path(p).name
        exp.file(p,'sync/20261022/controller_state_evidence/'+label,'shared',force=True)
    # A small actual input clip makes latest processing reproducible without publishing all recordings.
    samples=Path('/home/lzy/projects/external_connect/sample_videos')
    for p in samples.rglob('*'):
        if p.is_file():
            dest=Path('external_connect/sample_videos')/p.relative_to(samples)
            exp.file(p,dest,'samples',force=(p.name=='experiment_overlay_v1.mp4'))
    exp.file('/home/lzy/projects/microbubble_aggregation.md','docs/shared_project_overview_20261022.md','shared',force=True)
    exp.finish()
    print('elapsed_seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':main()
