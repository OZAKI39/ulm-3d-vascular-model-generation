"""Render existing facet WSS directly, using the current surface presentation."""
from pathlib import Path
import argparse
import json
import os
import sys
import time

os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW','vtkEGLRenderWindow')
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import numpy as np
import pyvista as pv
import imageio_ffmpeg
import render_surface_fields as base
import render_visualization as shared
from validate_results import validate_annotations

KIND='wss_raw'
base.SPECS[KIND]={**base.SPECS['wss'], 'scalar':'WSS_raw_Pa',
    'title':'Full vessel | Raw wall shear stress',
    'footer':'Raw WSS per wall triangle | No nodal averaging'}


class RawSurfaceView(base.SurfaceView):
    """Reuse all camera, geometry lighting, typography and overlay code."""
    def __init__(self,mesh,full_surface,orbit,ports,size=shared.SIZE,annotations=None):
        assert 'WSS_raw_Pa' in mesh.cell_data and 'WSS_raw_Pa' not in mesh.point_data
        super().__init__(mesh,full_surface,orbit,ports,KIND,size=size,annotations=annotations)
        mapper=self.actor.mapper
        # Select the named cell array explicitly; no cell-to-point filter.
        mapper.SetScalarModeToUseCellFieldData()
        mapper.SelectColorArray('WSS_raw_Pa')
        mapper.SetColorModeToMapScalars()
        mapper.InterpolateScalarsBeforeMappingOff()
        mapper.Update()
        assert mapper.GetScalarModeAsString()=='UseCellFieldData'
        assert mapper.GetArrayName()=='WSS_raw_Pa'
        dataset=pv.wrap(mapper.GetInput())
        assert dataset.n_cells==mesh.n_cells
        assert np.array_equal(dataset.cell_data['WSS_raw_Pa'],mesh.cell_data['WSS_raw_Pa'])
        assert np.array_equal(dataset.faces,mesh.faces)
        assert np.array_equal(dataset.points,mesh.points)
        self.mapper_record=dict(scalar_mode=mapper.GetScalarModeAsString(),
            scalar_name=mapper.GetArrayName(),interpolate_scalars_before_mapping=bool(mapper.GetInterpolateScalarsBeforeMapping()),
            cell_count=dataset.n_cells,cell_array_unchanged=True,geometry_and_connectivity_unchanged=True,
            normals_lighting='Original smooth normals and ambient/diffuse/specular preserved; no WSS averaging')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE)
    args=parser.parse_args();out=args.output.resolve();case=ROOT/'input_data'
    for name in ['figures','animations','inspection']:(out/name).mkdir(parents=True,exist_ok=True)
    movie=out/'animations/wss_raw_full_vessel_enlarged.mp4'
    if movie.exists():raise FileExistsError('Preserve completed output: '+str(movie))
    start=time.time();lock=json.loads((HERE/'audit/source_lock.json').read_text())
    for name,digest in lock.items():assert shared.sha(ROOT/name)==digest,name
    compute=json.loads((case/'field_diagnostics/COMPUTE_VALIDATION.json').read_text())
    assert compute['all_pass'] and compute['case']=='mean-2p0-mmps-A-ROI-only-balanced-pressure-v1'
    wall_path=case/'field_diagnostics/data/wall_wss_si.vtp'
    assert shared.sha(wall_path)==compute['outputs_sha256']['data/wall_wss_si.vtp']
    wall=pv.read(wall_path);raw=wall.cell_data['WSS_raw_Pa'].copy()
    areas=wall.cell_data['Area_m2'];traction=wall.cell_data['Tangential_viscous_traction_Pa']
    assert np.isfinite(raw).all() and np.all(areas>0)
    assert np.allclose(np.linalg.norm(traction,axis=1),raw,rtol=1e-13,atol=1e-13)
    assert len(raw)==compute['wall_facets'] and raw.min()>=0 and raw.max()<=55
    stats=dict(wall_triangles=wall.n_cells,raw_min_Pa=float(raw.min()),raw_max_Pa=float(raw.max()),
               area_weighted_mean_Pa=float(np.average(raw,weights=areas)),wall_area_m2=float(areas.sum()))
    assert np.isclose(stats['area_weighted_mean_Pa'],compute['wss_Pa']['area_weighted_mean'])
    wall.points=np.asarray(wall.points,dtype=float)*1e6
    full=pv.read(case/'field_diagnostics/data/pressure_surface_si.vtp')
    full.points=np.asarray(full.points,dtype=float)*1e6
    lines=pv.read(case/'streamlines/data/streamlines_si.vtp')
    framing=np.vstack([full.points,np.asarray(lines.points,dtype=float)*1e6])
    reference=json.loads((ROOT/'results/surface_fields/wss_camera.json').read_text())
    orbit=shared.Turntable(framing,[0,0,1],reference[0]['center_um'])
    assert orbit.trace==reference,'Camera must exactly match existing WSS movie'
    ports=[]
    for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        cap=pv.read(case/f'solver_mesh/mesh-surfaces/{role}.vtp')
        ports.append(('Inlet' if role=='INLET' else 'Outlet '+role[-2:],np.array(cap.center)*1e6))
    view=RawSurfaceView(wall,full,orbit,ports,size=(3840,2160))
    annotations=view.annotations
    for index in shared.KEYS[:4]:
        frame=view.frame(index,still=True)
        name='wss_raw_overview_4k' if index==0 else f'wss_raw_angle_{index*360/shared.FRAMES:03.0f}_4k'
        frame.save(out/'figures'/f'{name}.png',dpi=(300,300))
        if index==0:frame.convert('RGB').save(out/'figures'/f'{name}.pdf',resolution=300)
    annotations.export(out/'wss_raw_annotations.json')
    shared.dump(out/'wss_raw_camera.json',orbit.trace)
    assert (out/'wss_raw_annotations.json').read_bytes()==(ROOT/'results/surface_fields/wss_annotations.json').read_bytes()
    gpu=[s for s in view.p.render_window.ReportCapabilities().splitlines()
         if any(t in s for t in ['OpenGL vendor','OpenGL renderer','OpenGL version'])]
    assert any('NVIDIA GeForce RTX 4090' in s for s in gpu),gpu
    mapper_record=view.mapper_record;view.close()
    print('RAW_WSS_4K_STILLS_READY '+json.dumps(stats),flush=True)
    view=RawSurfaceView(wall,full,orbit,ports,annotations=annotations)
    assert view.mapper_record==mapper_record
    temporary=movie.with_suffix('.tmp.mp4')
    writer=imageio_ffmpeg.write_frames(str(temporary),size=shared.SIZE,fps=shared.FPS,
        codec='libx264',quality=8,pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=8,
        output_params=['-preset','medium','-movflags','+faststart','-threads','4'])
    writer.send(None)
    try:
        for i in range(shared.FRAMES):
            writer.send(np.ascontiguousarray(view.frame(i)))
            if i%72==0:print('raw_wss',i,'/',shared.FRAMES,flush=True)
    finally:
        writer.close();view.close()
    os.replace(temporary,movie)
    video=base.validate_video(movie,KIND,out)
    annotation_checks=validate_annotations(case,out,KIND,orbit.trace,
        surface_file='field_diagnostics/data/pressure_surface_si.vtp')
    for name,digest in lock.items():assert shared.sha(ROOT/name)==digest,name
    assert np.array_equal(wall.cell_data['WSS_raw_Pa'],raw)
    record=dict(all_pass=True,case=compute['case'],source_field_sha256=compute['source_field_sha256'],
        source_surface_sha256=shared.sha(wall_path),source_files_unchanged=True,source_lock=lock,
        raw_values_unchanged=True,source_fields_recomputed=False,cell_to_point_averaging=False,
        scalar_interpolation=False,scalar_clipping=False,geometry_smoothing=False,
        geometry_lighting_preserved=True,scalar='WSS_raw_Pa',scalar_association='CELL',
        colorbar_range_Pa=[0,55],stats=stats,mapper=mapper_record,video=video,
        camera_identical_to_existing_WSS=True,annotations_identical_to_existing_WSS=True,
        annotation_checks=annotation_checks,OpenGL=gpu,encoder='CPU libx264',
        new_CFD_calls=0,new_trajectory_integrations=0,script_sha256=shared.sha(Path(__file__)),
        elapsed_seconds=time.time()-start)
    shared.dump(out/'MEDIA_VALIDATION.json',record)
    print('RAW_WSS_RENDER_COMPLETE '+json.dumps(dict(all_pass=True,seconds=record['elapsed_seconds'])),flush=True)


if __name__=='__main__':main()
