"""Independently verify collected raw-facet animation and protected sources."""
from pathlib import Path
import hashlib,json
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import Image

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    record=json.loads((HERE/'MEDIA_VALIDATION.json').read_text())
    assert record['all_pass'] and record['scalar']=='WSS_raw_Pa' and record['scalar_association']=='CELL'
    assert not any(record[k] for k in ['source_fields_recomputed','cell_to_point_averaging',
        'scalar_interpolation','scalar_clipping','geometry_smoothing','new_CFD_calls','new_trajectory_integrations'])
    assert record['mapper']['scalar_mode']=='UseCellFieldData'
    assert record['mapper']['scalar_name']=='WSS_raw_Pa'
    assert not record['mapper']['interpolate_scalars_before_mapping']
    assert record['script_sha256']==sha(HERE/'render_raw_wss.py')
    for relative,digest in record['source_lock'].items():assert sha(ROOT/relative)==digest,relative
    assert (HERE/'wss_raw_camera.json').read_bytes()==(ROOT/'results/surface_fields/wss_camera.json').read_bytes()
    assert (HERE/'wss_raw_annotations.json').read_bytes()==(ROOT/'results/surface_fields/wss_annotations.json').read_bytes()
    surface=pv.read(ROOT/'input_data/field_diagnostics/data/wall_wss_si.vtp')
    raw=surface.cell_data['WSS_raw_Pa'];areas=surface.cell_data['Area_m2']
    assert np.allclose(np.linalg.norm(surface.cell_data['Tangential_viscous_traction_Pa'],axis=1),raw,rtol=1e-13,atol=1e-13)
    assert raw.min()==record['stats']['raw_min_Pa'] and raw.max()==record['stats']['raw_max_Pa']
    assert np.average(raw,weights=areas)==record['stats']['area_weighted_mean_Pa']
    assert surface.n_cells==record['mapper']['cell_count']==45221
    path=HERE/record['video']['file'];assert sha(path)==record['video']['sha256']
    reader=imageio_ffmpeg.read_frames(str(path),pix_fmt='rgb24')
    meta=next(reader);count=0;distinct=set()
    for frame in reader:count+=1;distinct.add(hashlib.sha256(frame).hexdigest())
    assert count==len(distinct)==432 and meta['fps']==24 and meta['duration']==18
    assert tuple(meta['size'])==(1920,1080)
    figures=[]
    for path in sorted((HERE/'figures').glob('*_4k.png')):
        with Image.open(path) as im:
            assert im.size==(3840,2160) and min(im.info['dpi'])>299
        figures.append(path.name)
    assert len(figures)==4
    result=dict(all_pass=True,source_files_verified=len(record['source_lock']),
        same_camera_and_labels_as_existing_WSS=True,original_movie_preserved=True,
        raw_cell_data_and_mapping_verified=True,raw_stats_Pa=record['stats'],
        video=dict(**meta,decoded_frames=count,distinct_frames=len(distinct)),figures_4k=figures,
        rendering_backend=record['OpenGL'],encoder=record['encoder'])
    (HERE/'LOCAL_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
