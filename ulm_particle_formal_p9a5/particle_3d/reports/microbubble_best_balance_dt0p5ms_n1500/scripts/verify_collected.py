"""Independent local check after incremental copy completes."""
from pathlib import Path
import hashlib,json,subprocess
import numpy as np
import imageio_ffmpeg
from PIL import Image
HERE=Path(__file__).resolve().parents[1]


def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def canonical_sha(v):
    return hashlib.sha256((json.dumps(v,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False)+'\n').encode()).hexdigest()


def main():
    delivery=json.loads((HERE/'DELIVERY_COMPLETE.json').read_text())
    for relative,expected in delivery['files'].items():assert sha(HERE/relative)==expected,relative
    cohort=json.loads((HERE/'data/cohort.json').read_text())
    summary=json.loads((HERE/'data/final_summary.json').read_text())
    native=json.loads((HERE/'data/native_exact_bytes_verification.json').read_text())
    assert native['PASS_FOR_ACTUAL_HANDOFF_DEPENDENCY']
    assert len(cohort['events'])==1500 and summary['dt_s']==.0005
    assert summary['integration_finished'] and summary['numerical_gate']
    assert summary['all_completion_hashes_verified'] and summary['GPU_diagnostics_pass']
    assert not any(summary['safety_totals'].values())
    assert summary['flow_sha256']=='fcc692caa74c70d7ecbf9ae6662d29d45abbae926b04783aa52886d382a4bdb2'
    identity=summary['identity']
    assert sha(HERE/'scripts/campaign.py')==identity['campaign_sha256']
    for name,digest in identity['acceleration_sources'].items():
        assert sha(HERE/'scripts'/name)==digest,name
    source=Path(json.loads((HERE/'config.json').read_text())['source_root'])
    if not source.exists():source=HERE.parent.parent.parent
    for relative,digest in identity['protected_sources'].items():
        assert sha(source/relative)==digest,relative
    count=samples=files=0;maxdt=0;statuses={};outlets={o:0 for o in ['O1','O2','O3']}
    for event in cohort['events']:
        folder=HERE/'tracks'/f"mb_{event['particle_id']:06d}"
        marker=json.loads((folder/'COMPLETE.json').read_text())
        assert marker['event_sha256']==canonical_sha(event)
        assert marker['identity']==summary['identity']
        for name,digest in marker['files'].items():
            path=(folder/name).resolve();assert path.is_relative_to(folder.resolve())
            assert sha(path)==digest,str(path);files+=1
        meta=json.loads((folder/'trajectory.json').read_text())
        metrics=json.loads((folder/'metrics.json').read_text())
        a=np.load(folder/'trajectory.npz')['samples']
        assert np.isfinite(a).all() and np.all(np.diff(a[:,0])>0)
        maxdt=max(maxdt,float(np.diff(a[:,0]).max()))
        assert meta['integration_config']['dt_s']==.0005
        assert np.array_equal(a[0,1:4],event['birth_center_m'])
        assert metrics['particle_id']==event['particle_id']
        statuses[metrics['status']]=statuses.get(metrics['status'],0)+1
        if metrics['outlet'] is not None:outlets[metrics['outlet']]+=1
        samples+=len(a);count+=1
    assert maxdt<=.0005+1e-12 and outlets==summary['outlets'] and statuses==summary['statuses']
    export=np.load(HERE/'data/trajectories_dt0p5ms.npz')
    assert export['offsets'].shape==(1501,) and np.all(np.diff(export['offsets'])>0)
    assert len(export['samples'])==export['offsets'][-1] and np.isfinite(export['samples']).all()
    images={}
    for path in sorted((HERE/'figures').glob('*.png')):
        with Image.open(path) as im:images[path.name]=dict(size=im.size,dpi=im.info.get('dpi'));assert min(im.size)>=1000
    decoder=imageio_ffmpeg.read_frames(str(HERE/'animations/microbubble_age_aligned.mp4'),pix_fmt='rgb24')
    video=next(decoder)
    video['decoded_frames']=sum(1 for frame in decoder)
    assert tuple(video['size'])==(1920,1080) and video['decoded_frames']==288 and video['fps']==24
    result=dict(PASS=True,tracks=count,raw_samples=samples,track_files_hash_verified=files,
        delivery_files_hash_verified=len(delivery['files']),max_accepted_dt_s=maxdt,
        protected_source_files_verified=len(identity['protected_sources']),
        campaign_and_acceleration_sources_verified=True,
        outlets=outlets,statuses=statuses,images=images,video=video,
        caveat='Integrity and explicit recorded numerical gates; not a time-convergence study')
    (HERE/'LOCAL_VERIFICATION.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=='__main__':main()
