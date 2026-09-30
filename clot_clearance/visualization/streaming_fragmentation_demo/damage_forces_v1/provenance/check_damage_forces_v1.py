"""Check saved-force visualization against original states and fully decode media."""
import argparse,hashlib,json,os
from pathlib import Path
import numpy as np
import imageio.v2 as imageio
import imageio_ffmpeg
from PIL import Image


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();out=a.output.resolve()
    m=json.loads((out/'RENDER_MANIFEST.json').read_text());z=np.load(Path(m['run'])/'states.npz');d=np.load(out/'force_samples.npz')
    checks={}
    checks['positions_unmodified']=np.array_equal(z['x'],d['positions_m'])
    checks['damage_unmodified']=np.array_equal(z['damage'],d['particle_damage'])
    checks['all_particles_in_every_state']=d['positions_m'].shape==(26,360,3)
    checks['external_force_sum']=np.array_equal(d['surface_force_N']+d['relaxation_force_N'],d['external_fluid_force_N'])
    checks['net_force_all_particles']=bool(np.allclose(d['external_fluid_force_N'].sum(axis=1),d['net_external_force_N'],rtol=1e-15,atol=1e-20))
    ids=d['displayed_arrow_particle_ids']
    checks['fixed_material_arrow_selection']=len(ids)==54 and len(np.unique(ids))==54 and not np.any(z['fixed'][ids])
    checks['single_global_linear_arrow_scale']=np.array_equal(d['displayed_arrow_vectors_mm'],d['external_fluid_force_N'][:,ids]*1e6*m['arrow_mm_per_uN'])
    checks['zero_force_on_base']=bool(np.all(d['external_fluid_force_N'][:,z['fixed']]==0))
    checks['read_only_force_provider_matches_saved_VTK']=m['checks']['maximum_surface_force_reconstruction_error_N']<1e-18 and m['checks']['maximum_saved_fluid_velocity_error_m_s']<1e-14
    checks['no_bubble_or_fragment_panel']=m['bubble_marker'] is False and m['fragment_ID_panel'] is False and m['bond_lines'] is False
    before=json.loads((out/'provenance/INPUT_FILES_SHA256.json').read_text())
    changed=[p for p,s in before.items() if not Path(p).is_file() or hashlib.sha256(Path(p).read_bytes()).hexdigest()!=s]
    checks['preexisting_files_unchanged']=not changed
    os.environ['IMAGEIO_FFMPEG_EXE']=imageio_ffmpeg.get_ffmpeg_exe()
    reader=imageio.get_reader(out/'particle_damage_external_forces.mp4');meta=reader.get_meta_data();count=0;size_ok=True;errors={}
    (out/'provenance/decoded_frames').mkdir(exist_ok=True)
    for i,frame in enumerate(reader):
        count+=1;size_ok &= frame.shape==(1080,1920,3)
        if i in [0,12,66,150]:
            k=i//6;ref=imageio.imread(out/f'frames/frame_{k:04d}.png')[:,:,:3]
            errors[str(k)]=float(np.abs(frame.astype(float)-ref).mean())
            imageio.imwrite(out/f'provenance/decoded_frames/frame_{k:04d}.png',frame)
    reader.close()
    checks['movie_full_decode']=count==156 and size_ok and meta['fps']==12 and max(errors.values())<3
    with Image.open(out/'particle_damage_external_forces.gif') as gif:
        n=gif.n_frames;duration=0
        for i in range(n):gif.seek(i);gif.load();duration+=gif.info.get('duration',0)
    checks['gif_full_decode']=n==26 and duration==13000
    result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,protected_files=len(before),changed_files=changed,
        video=dict(frames=count,size=meta['size'],fps=meta['fps'],duration_s=meta['duration'],event_frame_pixel_MAE=errors),
        gif=dict(frames=n,duration_ms=duration),
        force_units='N in NPZ; microNewton in animation and CSV',time_resolution='saved macro endpoints only; no subcycle reconstruction',
        physical_validation_claim=False)
    (out/'QC.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
