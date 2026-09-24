#!/usr/bin/env python3
"""Decode all exported bytes, verify frame-state receipts, then build storyboard."""
from pathlib import Path
import argparse,hashlib,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import imageio_ffmpeg
from PIL import Image,ImageDraw
from particle_3d.particle82_results import Scene
from particle_3d.particle82_visuals import ANIMATIONS
from particle_3d.particle81_replay import TITLE,FOOTERS
from particle_3d.particle8_full3d_visuals import text
from particle_3d.particle8_replay import read,digest,canonical_hash
from particle_3d.particle82_provenance import atomic_json


def audit(root):
    root=Path(root);scene=Scene(root);records=[]
    overview=Image.new('RGB',(1800,1850),'white');d=ImageDraw.Draw(overview)
    text(d,(30,15),'Particle-8.2 | 10 Storyboard from actual exported MP4s',28,bold=True)
    text(d,(30,56),TITLE,18,'#208ca6')
    for row,kind in enumerate(ANIMATIONS):
        name='particle8_2_anim_'+kind;manifest=read(root/'frames'/(name+'.json'));indices=manifest['keyframe_indices']
        video=root/'animations'/(name+'.mp4');reader=imageio_ffmpeg.read_frames(str(video),pix_fmt='rgb24');meta=next(reader)
        sheet=Image.new('RGB',(1800,450),'white');sd=ImageDraw.Draw(sheet);text(sd,(15,10),kind+' | actual decoded frames',23,bold=True)
        distinct=set();found=[];n=0
        for i,pixels in enumerate(reader):
            n=i+1;distinct.add(hashlib.sha256(pixels).digest());receipt=manifest['frames'][i]
            expected=scene.snapshot(receipt['physical_time_s'])
            if receipt['state']!=expected or canonical_hash(expected)!=receipt['state_sha256']:raise ValueError('Exported frame state drift')
            if i in indices:
                col=indices.index(i);im=Image.frombytes('RGB',meta['size'],pixels)
                im.save(root/'inspection'/f'{name}_decoded_{i:04d}.png');sheet.paste(im.resize((600,375)),(600*col,43))
                text(sd,(600*col+8,422),f'Frame {i} | t={receipt["physical_time_s"]:.6f} s',16)
                overview.paste(im.resize((580,363)),(col*600+10,row*400+95));found.append(i)
        gif_path=root/'animations'/(name+'.gif');gif=Image.open(gif_path)
        for k in range(gif.n_frames):gif.seek(k);gif.convert('RGB').load()
        sheet.save(root/'inspection'/(name+'_contact_sheet.png'))
        record=dict(kind=kind,mp4_sha256=digest(video),gif_sha256=digest(gif_path),decoded_mp4_frames=n,
            decoded_gif_frames=gif.n_frames,distinct_frames=len(distinct),size=meta['size'],fps=meta['fps'],codec=meta['codec'],
            sampled_indices=found,source_scene_sha256=scene.sha256,
            all_pass=n==manifest['frame_count'] and found==indices and len(distinct)>40 and gif.n_frames>30)
        records.append(record);print(kind,n,gif.n_frames,flush=True)
    for i,line in enumerate(FOOTERS):text(d,(30,1720+25*i),line,13,'#884b43')
    overview.save(root/'figures/10_particle8_2_storyboard.png')
    atomic_json(root/'figure_sources/10_particle8_2_storyboard.json',dict(source_scene_sha256=scene.sha256,
        method='FIRST_MIDDLE_LAST_FROM_ACTUAL_EXPORTED_VIDEO_BYTES',labels=[TITLE,*FOOTERS],media=records))
    atomic_json(root/'data/media_audit.json',dict(all_pass=all(r['all_pass'] for r in records),animations=records,
        method='ALL_MP4_AND_GIF_FRAMES_DECODED_AND_FRAME_RECEIPTS_RECOMPUTED',ffmpeg_version=imageio_ffmpeg.get_ffmpeg_version()))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);audit(p.parse_args().root)
