"""Hash deliverables and describe bulky remote assets without Git redistribution."""
import sys,json,gzip,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
def main():
    server='/workspace/f10w_vision_20261017';bundles=[]
    for file,remote,description in [('fullvideo_delivery.tar.gz',server+'/setup/fullvideo_delivery.tar.gz','All 937 frozen PIV chunk records and parquet outputs'),('models_complete.tar.gz',server+'/setup/models_complete.tar.gz','All model benchmark stages, real clips, two replay versions, solver diagnostics and evidence'),('logs_complete.tar.gz',server+'/setup/logs_complete.tar.gz','All F10W server logs and environment locks')]:
        p=OUT/'evidence'/file;bundles.append(dict(name=file,path=remote,bytes=p.stat().st_size,sha256=sha(p),description=description))
    files=[]
    for p in sorted((OUT/'results/full_video_chunks').glob('*')):
        if p.is_file():files.append(dict(local=str(p.relative_to(OUT)),remote=server+'/runs/full_video_e15b/'+str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=sha(p)))
    for folder in ['results/benchmark','results/replay_before_hash_cache_fix']:
        for p in sorted((OUT/folder).rglob('*')):
            if p.is_file():files.append(dict(local=str(p.relative_to(OUT)),remote=server+'/runs/models_v1/'+str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=sha(p)))
    for p in sorted((OUT/'results').glob('*.json')):
        if p.name.startswith('causal_replay_') or p.name=='real_comparison_clip.json':files.append(dict(local=str(p.relative_to(OUT)),remote=server+'/runs/models_v1/'+str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=sha(p)))
    dump(OUT/'REMOTE_ARTIFACTS.json',dict(server='ssh -p 15283 root@95.3.33.46',root=server,bundles=bundles,files=files,input_and_weight_manifest='evidence/remote_assets.json',old_failed_generator_directory=server+'/runs/models_v1/data_attempt_nonaff45',restore_script='scripts/restore_assets.py',source_snapshots=sorted(p.name for p in (OUT/'evidence/deployment').iterdir() if p.is_dir()),no_third_party_redistribution=True))
    manifest={}
    for folder in ['src','scripts','configs','tests','results','figures','videos','models','logs','evidence']:
        for p in sorted((OUT/folder).rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts or p.suffix in ['.tar','.gz'] or p.name.startswith('source.tar'):continue
            manifest[str(p.relative_to(OUT))]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for p in OUT.glob('*.md'):manifest[p.name]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for name in ['FINAL_GATE.json','FULL_VIDEO_PIV_SUMMARY.json','REMOTE_ARTIFACTS.json']:
        p=OUT/name;manifest[name]=dict(bytes=p.stat().st_size,sha256=sha(p))
    dump(OUT/'INPUT_OUTPUT_SHA256.json',dict(identity=identity(),input_source=dict(path=str(VIDEO),sha256=sha(VIDEO)),files=manifest,excluded='This manifest itself; Python environments, third-party working trees, transfer archives (separately hashed), unpublished large input caches indexed in remote_assets.json'))
    print('Manifest',len(manifest),'local files;',len(files),'remote bulk output references')
if __name__=='__main__':main()
