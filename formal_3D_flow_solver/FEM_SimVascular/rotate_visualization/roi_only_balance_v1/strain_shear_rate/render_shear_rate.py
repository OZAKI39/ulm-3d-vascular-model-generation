"""Current-style turntable of equivalent shear rate on the true exterior."""
from pathlib import Path
import argparse,json,os,sys,time
os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW','vtkEGLRenderWindow')
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import ImageFont
import render_surface_fields as base
import render_visualization as shared
from validate_results import validate_annotations

KIND='shear_rate';SCALAR='ShearRate_s_inv'

class RateView(base.SurfaceView):
    def __init__(self,mesh,full,orbit,ports,size=shared.SIZE,annotations=None):
        assert SCALAR in mesh.cell_data and SCALAR not in mesh.point_data
        super().__init__(mesh,full,orbit,ports,KIND,size=size,annotations=annotations)
        mapper=self.actor.mapper
        mapper.SetScalarModeToUseCellFieldData();mapper.SelectColorArray(SCALAR)
        mapper.SetColorModeToMapScalars();mapper.InterpolateScalarsBeforeMappingOff();mapper.Update()
        actual=pv.wrap(mapper.GetInput())
        assert mapper.GetScalarModeAsString()=='UseCellFieldData' and mapper.GetArrayName()==SCALAR
        assert np.array_equal(actual.cell_data[SCALAR],mesh.cell_data[SCALAR])
        assert np.array_equal(actual.points,mesh.points) and np.array_equal(actual.faces,mesh.faces)
        self.mapping=dict(scalar=mapper.GetArrayName(),mode=mapper.GetScalarModeAsString(),
            interpolate_scalars_before_mapping=bool(mapper.GetInterpolateScalarsBeforeMapping()),
            cell_values_unchanged=True,topology_and_coordinates_unchanged=True,cells=actual.n_cells)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE)
    args=parser.parse_args();out=args.output.resolve();case=ROOT/'input_data';start=time.time()
    for d in ['figures','animations','inspection']:(out/d).mkdir(parents=True,exist_ok=True)
    movie=out/'animations/shear_rate_full_vessel_enlarged.mp4'
    if movie.exists():raise FileExistsError('Preserve completed MP4: '+str(movie))
    compute=json.loads((HERE/'COMPUTE_VALIDATION.json').read_text());assert compute['all_pass']
    for name,digest in compute['source_lock'].items():assert shared.sha(ROOT/name)==digest,name
    for name,digest in compute['files'].items():assert shared.sha(HERE/name)==digest,name
    base.SPECS[KIND]=dict(title='Full vessel | Shear rate magnitude',scalar=SCALAR,
        unit='Rate (s⁻¹)',clim=tuple(compute['colorbar_range_s_inv']),
        ticks=[int(x) for x in compute['colorbar_ticks_s_inv']],
        footer='Equivalent shear rate | Exterior facets from adjacent FEM cells',file='shear_rate_exterior_si.vtp')
    font=ImageFont.truetype(shared.FONT,round(shared.SIZE[1]*.021))
    assert .904*1920+font.getlength(base.SPECS[KIND]['unit'])<1910
    tick_x=round(.922*1920)+round(.014*1920)+.010*1920
    assert tick_x+font.getlength(str(int(compute['colorbar_range_s_inv'][1])))<1910
    mesh=pv.read(HERE/'data/shear_rate_exterior_si.vtp');raw=mesh.cell_data[SCALAR].copy()
    assert np.isfinite(raw).all() and raw.min()>=0 and raw.max()<=base.SPECS[KIND]['clim'][1]
    mesh.points=np.asarray(mesh.points,dtype=float)*1e6
    full=pv.read(case/'field_diagnostics/data/pressure_surface_si.vtp');full.points=np.asarray(full.points,dtype=float)*1e6
    assert np.array_equal(mesh.points,full.points) and np.array_equal(mesh.faces,full.faces)
    lines=pv.read(case/'streamlines/data/streamlines_si.vtp')
    framing=np.vstack([full.points,np.asarray(lines.points,dtype=float)*1e6])
    previous=json.loads((ROOT/'results/surface_fields/wss_camera.json').read_text())
    orbit=shared.Turntable(framing,[0,0,1],previous[0]['center_um']);assert orbit.trace==previous
    ports=[]
    for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        cap=pv.read(case/f'solver_mesh/mesh-surfaces/{role}.vtp')
        ports.append(('Inlet' if role=='INLET' else 'Outlet '+role[-2:],np.array(cap.center)*1e6))
    view=RateView(mesh,full,orbit,ports,size=(3840,2160));annotations=view.annotations
    for i in shared.KEYS[:4]:
        frame=view.frame(i,still=True)
        stem='shear_rate_overview_4k' if i==0 else f'shear_rate_angle_{i*360/shared.FRAMES:03.0f}_4k'
        frame.save(out/'figures'/f'{stem}.png',dpi=(300,300))
        if i==0:frame.convert('RGB').save(out/'figures'/f'{stem}.pdf',resolution=300)
    annotations.export(out/'shear_rate_annotations.json');shared.dump(out/'shear_rate_camera.json',orbit.trace)
    assert (out/'shear_rate_annotations.json').read_bytes()==(ROOT/'results/surface_fields/wss_annotations.json').read_bytes()
    gpu=[s for s in view.p.render_window.ReportCapabilities().splitlines()
         if any(t in s for t in ['OpenGL vendor','OpenGL renderer','OpenGL version'])]
    assert any('NVIDIA GeForce RTX 4090' in s for s in gpu),gpu
    mapping=view.mapping;view.close();print('SHEAR_RATE_4K_READY',flush=True)
    view=RateView(mesh,full,orbit,ports,annotations=annotations);assert view.mapping==mapping
    tmp=movie.with_suffix('.tmp.mp4')
    writer=imageio_ffmpeg.write_frames(str(tmp),size=shared.SIZE,fps=shared.FPS,codec='libx264',quality=8,
        pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=8,
        output_params=['-preset','medium','-movflags','+faststart','-threads','4'])
    writer.send(None)
    try:
        for i in range(shared.FRAMES):
            writer.send(np.ascontiguousarray(view.frame(i)))
            if i%72==0:print('shear_rate',i,'/',shared.FRAMES,flush=True)
    finally:writer.close();view.close()
    os.replace(tmp,movie);video=base.validate_video(movie,KIND,out)
    checks=validate_annotations(case,out,KIND,orbit.trace,'field_diagnostics/data/pressure_surface_si.vtp')
    for name,digest in compute['source_lock'].items():assert shared.sha(ROOT/name)==digest,name
    for name,digest in compute['files'].items():assert shared.sha(HERE/name)==digest,name
    assert np.array_equal(mesh.cell_data[SCALAR],raw)
    result=dict(all_pass=True,case=compute['source_case'],source_field_sha256=compute['source_field_sha256'],
        scalar=SCALAR,unit='s^-1',definition=compute['definition'],surface='Entire exterior including inlet/outlet caps',
        rendered_range_s_inv=[float(raw.min()),float(raw.max())],colorbar_range_s_inv=compute['colorbar_range_s_inv'],
        mapping=mapping,cell_to_point_averaging=False,scalar_clipping=False,geometry_smoothing=False,
        original_smooth_normal_lighting_preserved=True,camera_identical_to_existing_WSS=True,
        annotations_identical_to_existing_WSS=True,annotation_checks=checks,
        source_files_unchanged=True,derived_data_unchanged=True,video=video,OpenGL=gpu,encoder='CPU libx264',
        CFD_calls=0,new_trajectory_integrations=0,compute_record_sha256=shared.sha(HERE/'COMPUTE_VALIDATION.json'),
        script_sha256=shared.sha(Path(__file__)),elapsed_seconds=time.time()-start)
    shared.dump(out/'MEDIA_VALIDATION.json',result)
    print('SHEAR_RATE_RENDER_COMPLETE '+json.dumps(dict(all_pass=True,seconds=result['elapsed_seconds'])),flush=True)

if __name__=='__main__':main()
