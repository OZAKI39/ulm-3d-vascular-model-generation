#!/usr/bin/env python3
"""Decode every exported frame; construct Figure 08 from actual MP4 bytes."""
from pathlib import Path
import hashlib,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import imageio_ffmpeg
from PIL import Image,ImageDraw
from particle_3d.particle81_simulation import OUTPUT,dump
from particle_3d.particle81_replay import Scene,TITLE,FOOTERS
from particle_3d.particle81_visuals import ANIMATIONS
from particle_3d.particle8_full3d_visuals import text
from particle_3d.particle8_replay import read,digest


def main():
    root=OUTPUT;scene=Scene();rows=[]
    overview=Image.new('RGB',(1800,1850),'white');d=ImageDraw.Draw(overview)
    text(d,(30,15),'Particle-8.1 | 08 Exported-animation storyboard',28,bold=True)
    text(d,(30,56),TITLE,18,'#208ca6')
    for row,kind in enumerate(ANIMATIONS):
        name='particle8_1_anim_'+kind;manifest=read(root/'frames'/(name+'_frames.json'));indices=manifest['keyframe_indices']
        reader=imageio_ffmpeg.read_frames(str(root/'animations'/(name+'.mp4')),pix_fmt='rgb24');meta=next(reader)
        sheet=Image.new('RGB',(1800,450),'white');sd=ImageDraw.Draw(sheet);text(sd,(15,10),kind+' | actual decoded MP4 frames',23,bold=True)
        distinct=set();found=[]
        for i,pixels in enumerate(reader):
            distinct.add(hashlib.sha256(pixels).digest())
            if i in indices:
                col=indices.index(i);im=Image.frombytes('RGB',meta['size'],pixels)
                im.save(root/'inspection'/f'{name}_decoded_{i:04d}.png');sheet.paste(im.resize((600,375)),(col*600,43))
                f=manifest['frames'][i];text(sd,(col*600+8,422),f'Frame {i} | t={f["physical_time_s"]:.6f} s',16)
                overview.paste(im.resize((580,363)),(col*600+10,row*400+95));found.append(i)
        gif=Image.open(root/'animations'/(name+'.gif'))
        for k in range(gif.n_frames):gif.seek(k);gif.convert('RGB').load()
        sheet.save(root/'inspection'/(name+'_contact_sheet.png'))
        rows.append(dict(kind=kind,mp4_sha256=digest(root/'animations'/(name+'.mp4')),gif_sha256=digest(root/'animations'/(name+'.gif')),
            decoded_mp4_frames=i+1,decoded_gif_frames=gif.n_frames,distinct_frames=len(distinct),size=meta['size'],fps=meta['fps'],
            codec=meta['codec'],sampled_indices=found,source_scene_sha256=scene.sha256,
            pass_checks=(i+1)==manifest['frame_count'] and found==indices and len(distinct)>40 and gif.n_frames>30))
        print(kind,i+1,gif.n_frames,flush=True)
    for i,line in enumerate(FOOTERS):text(d,(30,1720+25*i),line,13,'#884b43')
    overview.save(root/'figures/particle8_1_fig_08_exported_animation_storyboard.png')
    dump(root/'data/particle8_1_fig_08_exported_animation_storyboard_source.json',dict(source_scene_sha256=scene.sha256,
        figure=8,details='FIRST_MIDDLE_LAST_DECODED_FROM_ALL_ACTUAL_MP4_FILES',scientific_labels=[TITLE,*FOOTERS],media=rows))
    dump(root/'data/media_audit.json',dict(all_pass=all(r['pass_checks'] for r in rows),animations=rows,method='ALL_MP4_AND_GIF_FRAMES_DECODED',
        ffmpeg_version=imageio_ffmpeg.get_ffmpeg_version()))


if __name__=='__main__':main()
