"""V4 audit, copied from v3: component identity, interpolation and decoded video."""
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
        no_center_marker_or_force_chart=not any(m[k] for k in ['spatial_center_marker','total_particle_force_arrows','force_chart','force_numbers','cycle_label','bubble_sphere']),
        Arial_regular=m['font']['family']=='Arial' and m['font']['weight']=='normal' and sha(Path(m['font']['path']))==m['font']['sha256'])
    c=json.loads((Path(m['run'])/'CONFIG.json').read_text())
    t=np.load(out/'traction_samples.npz');raw_t=t['local_streaming_traction_Pa'];normal=t['surface_normals'];mask=t['load_eligible']
    expected_mask=z['surface'] & z['attached'] & ~z['fixed'][None,:] & (np.linalg.norm(normal,axis=2)>0)
    checks['eligibility_matches_actual_surface_loading']=np.array_equal(mask,expected_mask)
    s=c['streaming'];points=z['x'];times=t['mechanical_time_s']
    envelope=np.exp(-.5*np.sum((points-np.array(s['bubble_center_m']))**2,axis=2)/s['streaming_length_scale_m']**2)
    modulation=s['mean_fraction']+s['oscillatory_amplitude']*np.sin(2*np.pi*c['simulation']['representative_frequency_Hz']*times+s['phase_rad'])
    tangent=np.array([1.,0.,0.])-normal[:,:,0,None]*normal
    formula=envelope[:,:,None]*modulation[:,None,None]*s['traction_scale_Pa']*tangent
    formula[~mask]=0
    checks['independent_local_traction_formula']=bool(np.allclose(raw_t,formula,rtol=1e-13,atol=1e-12))
    component_error=float(np.max(np.abs(raw_t+t['background_pipe_traction_Pa']-t['saved_total_surface_traction_Pa'])))
    normal_error=float(np.max(np.abs(np.sum(raw_t*normal,axis=2))))
    checks['local_plus_pipe_reproduces_saved_total']=component_error<1e-10
    checks['original_local_traction_is_tangential']=normal_error<1e-10
    checks['no_traction_on_unloaded_or_free_particles']=bool(np.all(raw_t[~mask]==0))
    checks['traction_exact_at_original_keyframes']=np.array_equal(d['display_streaming_traction_Pa'][key],raw_t)
    common=mask[left]&mask[right]
    expected_t=np.where(common[:,:,None],(1-alpha[:,None,None])*raw_t[left]+alpha[:,None,None]*raw_t[right],raw_t[left]);expected_t[key]=raw_t
    checks['traction_uses_documented_display_interpolation']=np.array_equal(d['display_streaming_traction_Pa'],expected_t)
    ids=t['displayed_arrow_particle_ids']
    checks['fixed_material_ids_and_masks']=np.array_equal(d['displayed_arrow_particle_ids'],ids) and np.array_equal(d['arrow_visible'],mask[left][:,ids])
    checks['one_global_linear_arrow_scale']=np.array_equal(d['arrow_vectors_mm'],expected_t[:,ids]*m['traction_scale_mm_per_Pa'])
    expected_normals=(1-alpha[:,None,None])*normal[left]+alpha[:,None,None]*normal[right]
    expected_normals/=np.maximum(np.linalg.norm(expected_normals,axis=2,keepdims=True),1e-30)
    expected_normals=np.where(common[:,:,None],expected_normals,normal[left]);expected_normals[key]=normal
    expected_tails=d['display_positions_m'][:,ids]*1000+m['glyph_outward_offset_mm']*expected_normals[:,ids]
    checks['glyph_offset_is_only_for_display']=np.array_equal(d['arrow_tails_mm'],expected_tails)
    before=json.loads((out/'provenance/INPUT_FILES_SHA256.json').read_text());changed=[p for p,h in before.items() if not Path(p).is_file() or sha(Path(p))!=h]
    checks['all_original_files_unchanged']=not changed
    step=np.linalg.norm(np.diff(d['display_positions_m'],axis=0),axis=2).max(axis=1)
    raw=np.linalg.norm(np.diff(z['x'],axis=0),axis=2).max(axis=1)
    checks['maximum_display_step_reduced_30x']=bool(np.isclose(step.max(),raw.max()/30,rtol=1e-10,atol=1e-15))
    ffmpeg=imageio_ffmpeg.get_ffmpeg_exe();os.environ['IMAGEIO_FFMPEG_EXE']=ffmpeg
    movie=out/'particle_damage_streaming_traction_60fps.mp4';reader=imageio.get_reader(movie);meta=reader.get_meta_data()
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
    with Image.open(out/'particle_damage_streaming_traction_20fps.gif') as gif:
        frames=gif.n_frames;duration=0
        for i in range(frames):gif.seek(i);gif.load();duration+=gif.info.get('duration',0)
    checks['GIF_full_decode']=duration==12550 and 1<frames<=251
    result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,protected_files=len(before),changed_files=changed,
        video=dict(frame_count=count,fps=meta['fps'],duration_s=meta['duration'],size=meta['size'],identical_adjacent_decoded_frames=duplicates,
            keyframe_pixel_MAE=MAE,minimum_timestamp_step_s=float(np.diff(pts).min()),maximum_timestamp_step_s=float(np.diff(pts).max())),
        traction=dict(maximum_decomposition_error_Pa=component_error,maximum_normal_component_Pa=normal_error,peak_saved_traction_Pa=float(np.linalg.norm(raw_t,axis=2).max()),sampled_material_ids=len(ids),visible_arrow_counts_at_saved_states=d['arrow_visible'][key].sum(axis=1).tolist()),
        motion=dict(maximum_original_snapshot_step_mm=float(raw.max()*1000),maximum_display_frame_step_mm=float(step.max()*1000)),
        gif=dict(frames=frames,duration_ms=duration),
        note='60 fps is display-only. Saved local surface traction excludes pipe stress and fragment relaxation. Vector interpolation does not reconstruct unrecorded subcycle oscillation; damage and loading masks retain macro-step changes.')
    (out/'QC.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
