"""Audit interpolation, preservation, actual decoded frames and constant cadence."""
import argparse,hashlib,json,os,re,subprocess
from pathlib import Path
import numpy as np
import imageio.v2 as imageio
import imageio_ffmpeg
from PIL import Image


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();out=a.output.resolve()
    m=json.loads((out/'RENDER_MANIFEST.json').read_text());z=np.load(Path(m['run'])/'states.npz');d=np.load(out/'display_samples.npz')
    key=np.array(m['original_endpoint_indices']);left=d['source_state_left'];right=d['source_state_right'];alpha=d['interpolation_alpha']
    expected=(1-alpha[:,None,None])*z['x'][left]+alpha[:,None,None]*z['x'][right];expected[:,z['fixed']]=z['X'][z['fixed']];expected[key]=z['x']
    checks=dict(original_positions_exact_at_keyframes=np.array_equal(d['display_positions_m'][key],z['x']),
        original_damage_exact_at_keyframes=np.array_equal(d['display_damage'][key],z['damage']),
        display_positions_are_only_linear_interpolation=np.array_equal(d['display_positions_m'],expected),
        damage_never_interpolated=np.array_equal(d['display_damage'],z['damage'][left]),
        anchors_remain_exact=bool(np.all(d['display_positions_m'][:,z['fixed']]==z['X'][z['fixed']])),
        all_particles_preserved=d['display_positions_m'].shape==(751,360,3),
        no_particle_arrow_or_force_chart=not any(m[k] for k in ['force_arrows','force_chart','force_numbers','cycle_label','bubble_sphere']),
        Arial_regular=m['font']['family']=='Arial' and m['font']['weight']=='normal' and sha(Path(m['font']['path']))==m['font']['sha256'])
    c=json.loads((Path(m['run'])/'CONFIG.json').read_text())
    checks['load_source_from_original_configuration']=np.array_equal(d['source_position_m'],c['streaming']['bubble_center_m'])
    checks['source_upstream_of_clot']=d['source_position_m'][0]<c['clot']['origin_m'][0]
    checks['explicit_actual_damage_color_range']=m['damage_color_range'][0]==0 and m['damage_color_range'][1]>=float(z['damage'].max()) and m['damage_color_range'][1]<.001
    before=json.loads((out/'provenance/INPUT_FILES_SHA256.json').read_text());changed=[p for p,h in before.items() if not Path(p).is_file() or sha(Path(p))!=h]
    checks['all_original_files_unchanged']=not changed
    step=np.linalg.norm(np.diff(d['display_positions_m'],axis=0),axis=2).max(axis=1)
    raw=np.linalg.norm(np.diff(z['x'],axis=0),axis=2).max(axis=1)
    checks['maximum_display_step_reduced_30x']=bool(np.isclose(step.max(),raw.max()/30,rtol=1e-10,atol=1e-15))
    ffmpeg=imageio_ffmpeg.get_ffmpeg_exe();os.environ['IMAGEIO_FFMPEG_EXE']=ffmpeg
    movie=out/'regularized_damage_transport.mp4';reader=imageio.get_reader(movie);meta=reader.get_meta_data()
    folder=out/'provenance/decoded_frames';folder.mkdir(exist_ok=True)
    count=0;shape=True;previous=None;duplicates=0;MAE={};black=True
    for i,frame in enumerate(reader):
        count+=1;shape &= frame.shape==(1080,1920,3);black &= np.max(frame[:10,:10])<=1
        if previous is not None:duplicates+=int(np.array_equal(frame,previous))
        previous=frame.copy()
        if i in [0,60,330,750]:
            source=imageio.imread(out/f'keyframes/state_{i//30:04d}.png')[:,:,:3]
            MAE[str(i)]=float(np.abs(source.astype(float)-frame).mean());imageio.imwrite(folder/f'frame_{i:04d}.png',frame)
    reader.close()
    checks['complete_60fps_video']=count==751 and meta['fps']==60 and shape and max(MAE.values())<3
    checks['pure_black_background']=bool(black)
    probe=subprocess.run([ffmpeg,'-hide_banner','-i',str(movie),'-vf','showinfo','-an','-f','null','-'],capture_output=True,text=True,check=True)
    pts=np.array([float(v) for v in re.findall(r'\bpts_time:([0-9.eE+-]+)',probe.stderr)])
    checks['constant_frame_timestamps']=len(pts)==751 and bool(np.allclose(np.diff(pts),1/60,rtol=0,atol=7e-5))
    with Image.open(out/'regularized_damage_transport.gif') as gif:
        frames=gif.n_frames;duration=0
        for i in range(frames):gif.seek(i);gif.load();duration+=gif.info.get('duration',0)
    checks['GIF_full_decode']=duration==12550 and 1<frames<=251
    checks={key:bool(value) for key,value in checks.items()}
    result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,protected_files=len(before),changed_files=changed,
        video=dict(frame_count=count,fps=meta['fps'],duration_s=meta['duration'],size=meta['size'],identical_adjacent_decoded_frames=duplicates,
            keyframe_pixel_MAE=MAE,minimum_timestamp_step_s=float(np.diff(pts).min()),maximum_timestamp_step_s=float(np.diff(pts).max())),
        motion=dict(maximum_original_snapshot_step_mm=float(raw.max()*1000),maximum_display_frame_step_mm=float(step.max()*1000)),
        gif=dict(frames=frames,duration_ms=duration),
        note='60 fps is a display product. Interpolated positions are not additional solver samples; damage retains original macro-step changes.')
    (out/'QC.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
