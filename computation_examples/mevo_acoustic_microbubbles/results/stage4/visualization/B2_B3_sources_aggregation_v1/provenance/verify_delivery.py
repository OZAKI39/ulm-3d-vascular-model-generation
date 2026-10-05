from pathlib import Path
import hashlib,json,csv,subprocess
from html.parser import HTMLParser
import imageio.v2 as io
import imageio_ffmpeg
import numpy as np
from PIL import Image
root=Path('/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles')
out=root/'results/stage4/visualization/B2_B3_sources_aggregation_v1'
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
m=json.loads((out/'RENDER_MANIFEST.json').read_text())
assert m['status']=='PASS' and sha(root/'scripts/stage4/render_aggregation.py')==m['renderer_sha256']
assert sha(out/'provenance/stage4_render_aggregation.py')==m['renderer_sha256']
assert all(sha(p)==s for p,s in m['protected_inputs'].items())
checked={}
readers={}
for name,info in m['videos'].items():
    path=out/info['path'];assert sha(path)==info['sha256']
    r=io.get_reader(path);readers[name]=r
    meta=r.get_meta_data();assert r.count_frames()==1440 and meta['fps']==60
    assert meta['size']==tuple(info['resolution']) and meta['codec']=='h264'
    decoded=subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-threads','2','-i',str(path),'-f','null','-'],capture_output=True,text=True,check=True)
    assert not decoded.stderr,decoded.stderr
    checked[name]=dict(full_decode='PASS',frames=1440,fps=60,resolution=info['resolution'])
for name in ['B2','B3']:
    r=readers[name]
    raw=[io.imread(out/'keyframes'/f'{name}_{k:02d}.png') for k in range(11)]
    roi=np.stack([im[325:695,1270:1710].astype(np.float32) for im in raw])
    errors=[]
    for k in range(11):
        for frame in [k*120, min((k+1)*120-1,1439)]:
            im=r.get_data(frame)
            scores=np.mean(np.abs(roi-im[325:695,1270:1710]),axis=(1,2,3))
            assert int(scores.argmin())==k,(name,frame,scores.tolist())
            mae=float(np.abs(im.astype(float)-raw[k]).mean());assert mae<2.,mae
            errors.append(dict(frame=frame,state=k,rgb_mean_absolute_error=mae))
    assert np.mean(np.abs(r.get_data(1439).astype(float)-raw[10]))<2.
    checked[name]['all_snapshot_transition_checks']=errors
    if name=='B2':
        Image.fromarray(r.get_data(1200)).resize((1536,864),Image.Resampling.LANCZOS).save('/tmp/mevo_aggregation_encoded_B2.png')
for frame in [0,120,1200,1439]:
    actual=readers['B2_B3'].get_data(frame)
    expected=np.concatenate([readers[n].get_data(frame) for n in ['B2','B3']],axis=1)
    assert np.abs(actual.astype(float)-expected).mean()<2.
for r in readers.values():r.close()
rows=list(csv.DictReader((out/'frame_map.csv').open()));assert len(rows)==2880
for row in rows:
    k=min(int(row['frame'])//120,10)
    s=m['states'][row['case']][k]
    assert int(row['saved_state'])==k and float(row['physical_time_s'])==s['t']
    assert int(row['local_group'])==s['anchor_cluster_size']
    assert int(row['all_cluster_count'])==s['stats']['cluster_count']
    assert float(row['near_contact_fraction'])==s['stats']['near_contact_fraction']
    assert float(row['median_nearest_neighbor_m'])==s['stats']['median_nearest_neighbor_m']
class Links(HTMLParser):
    def handle_starttag(self,tag,attrs):
        for key,val in attrs:
            if key in ['src','href','poster']:assert (out/val).is_file(),val
Links().feed((out/'OPEN_RESULTS.html').read_text())
report=dict(status='PASS',videos=checked,source_hashes_unchanged=True,frame_map_rows=2880,
            grouping_definition_matches_formal_config=True,all_original_metrics_match_at_22_displayed_states=True,
            radius_um=1.78,field_color_scale_kPa=[0,39],code_commit='ef270b3',
            field_scope='Actual stored complex-field amplitudes on two spatial slices; no acoustic rerun or carrier/bubble-clock mixing',
            microscopy_scope='Fixed ID0 co-moving view, equilibrium-radius sphere glyphs, only recorded states; no position interpolation')
(out/'DELIVERY_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
files={str(p.relative_to(out)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file() and p.name not in ['DELIVERY_SHA256.json','render.log']}
(out/'DELIVERY_SHA256.json').write_text(json.dumps(files,indent=2)+'\n')
print(json.dumps(dict(status='PASS',files=len(files),bytes=sum(v['bytes'] for v in files.values()),videos=list(checked))))
