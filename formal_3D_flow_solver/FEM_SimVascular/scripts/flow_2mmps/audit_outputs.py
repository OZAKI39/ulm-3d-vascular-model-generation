"""Decode actual MP4 files and independently audit display/frozen artifacts."""
from pathlib import Path
import json, subprocess, re, hashlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import imageio_ffmpeg
from validate import CASE,ROOT,sha,dump

def image_check(path):
    with Image.open(path) as im:
        im.load();array=np.asarray(im.convert('RGB'));width,height=im.size
    assert width>=1200 and height>=650, 'Figure too small'
    assert np.max(array[:8,:8])<=3 and np.max(array[-8:,-8:])<=3, 'Background is not black'
    foreground=np.max(array,axis=2)>25
    assert foreground.mean()>.004, 'Blank or nearly empty image'
    return dict(status='PASS',file=str(path.relative_to(CASE)) if path.is_relative_to(CASE) else str(path),
                sha256=sha(path),width=width,height=height,foreground_fraction=float(foreground.mean()))

def decode_video(path,expected_frames=432):
    ffmpeg=imageio_ffmpeg.get_ffmpeg_exe()
    cmd=[ffmpeg,'-v','error','-xerror','-i',str(path),'-map','0:v:0','-progress','pipe:1','-nostats','-f','null','-']
    result=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
    assert result.returncode==0 and not result.stderr.strip(), 'Animation decode failed: '+result.stderr[-1000:]
    counts=re.findall(r'^frame=(\d+)',result.stdout,re.M)
    assert counts and int(counts[-1])==expected_frames, 'Wrong decoded frame count'
    reader=imageio_ffmpeg.read_frames(str(path));meta=next(reader)
    assert tuple(meta['size'])==(1920,1080) and abs(meta['fps']-24)<1e-9
    assert abs(meta['duration']-18)<.05
    selected={};frames=[0,144,288]
    for index,data in enumerate(reader):
        if index in frames:
            a=np.frombuffer(data,np.uint8).reshape(1080,1920,3).copy()
            assert np.max(a[:8,:8])<12 and np.max(a[-8:,-8:])<12
            body=a[100:950,100:1570].astype(float)
            assert (body.max(2)-body.min(2)>20).sum()>1000, 'Missing colored vessel'
            right=a[200:850,1650:1850].astype(float)
            assert (right.max(2)-right.min(2)>30).sum()>3000, 'Missing right colorbar'
            selected[index]=a
        if index==frames[-1]:break
    reader.close()
    assert len(selected)==3
    assert np.mean(np.abs(selected[0].astype(float)-selected[144].astype(float)))>.5, 'No visible camera rotation'
    return dict(status='PASS',file=str(path.relative_to(CASE)) if path.is_relative_to(CASE) else str(path),
        sha256=sha(path),metadata=meta,decoded_frames=int(counts[-1]),decoder_exit_code=result.returncode,
        decoder_stderr=result.stderr,decode_command=cmd,selected_frames=frames),selected

def main():
    reports=CASE/'reports';manifest=json.loads((reports/'render_manifest.json').read_text())
    assert manifest['source_flow_sha256']==sha(CASE/'frozen_flow/steady_flow_mean_2p0_mmps.vtu')
    videos=[];panels=[]
    for i,item in enumerate(manifest['animations']):
        path=CASE/'animations'/item['file'];assert sha(path)==item['sha256']
        record,selected=decode_video(path);videos.append(record)
        trace=json.loads(path.with_name(path.stem+'_camera.json').read_text())
        assert len(trace)==432 and np.allclose(np.diff([t['azimuth_deg'] for t in trace]),360/432)
        assert np.allclose(np.diff([t['time_s'] for t in trace]),1/24)
        assert len({t['parallel_scale_um'] for t in trace})==1
        if i==0:assert max(abs(v) for t in trace for v in t['projected_bounds'])<.94
        for index,a in selected.items():
            panels.append((i,index,Image.fromarray(a).resize((960,540),Image.Resampling.LANCZOS)))
    assert len(videos)==2
    board=Image.new('RGB',(2880,1280),'black');draw=ImageDraw.Draw(board)
    font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',30)
    draw.text((30,20),'Rotation storyboard | decoded MP4 frames',fill='white',font=font)
    for j,(video,index,im) in enumerate(panels):
        x=(j%3)*960;y=110+(j//3)*610
        draw.text((x+28,y-43),f"{'Full vessel' if video==0 else 'Branch detail'} | {index/24:.1f} s | {45+index*360/432:.0f} deg",font=font,fill='#c9d3e0')
        board.paste(im,(x,y))
    storyboard=CASE/'figures/Figure_04_storyboard.png';board.save(storyboard)
    figures=[image_check(p) for p in sorted((CASE/'figures').glob('Figure_*.png'))]
    assert len(figures)>=6
    ffmpeg=imageio_ffmpeg.get_ffmpeg_exe()
    record=dict(status='PASS',figures=figures,animations=videos,
        storyboard=dict(path=str(storyboard.relative_to(CASE)),source='actual decoded MP4 frames',frames=[0,144,288]),
        full_vessel_all_angles_inside_viewport=True,camera_uniform_angular_increment=True,
        ffmpeg_version=subprocess.check_output([ffmpeg,'-version'],text=True).splitlines()[0],ffmpeg_sha256=sha(ffmpeg))
    dump(reports/'artifact_validation.json',record)
    print(json.dumps(dict(status='PASS',figures=len(figures),animations=len(videos),decoded_frames=[v['decoded_frames'] for v in videos]),indent=2))

if __name__=='__main__':main()
