"""Independent data, camera, label, and presentation checks for pressure/WSS."""
import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pyvista as pv

from validate_results import sha, validate_annotations

HERE = Path(__file__).resolve().parent


def validate(case,output):
    media = json.loads((output/'MEDIA_VALIDATION.json').read_text())
    assert media['all_pass'] and media['source_values_unchanged']
    assert media['script_sha256'] == sha(HERE/'render_surface_fields.py')
    assert media['shared_renderer_sha256'] == sha(HERE/'render_visualization.py')
    for file in ['render_surface_fields.py','render_visualization.py']:
        tree = ast.parse((HERE/file).read_text())
        assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
                       and n.func.attr=='rectangle' for n in ast.walk(tree))
    for name,value in media['source_sha256'].items():
        assert sha(case/name) == value
    expected = {
        'pressure': ('Steady FEM pressure on the vessel surface','Pressure (Pa)',[-5.,5000.]),
        'wss': ('Derived from FEM velocity gradient','WSS (Pa)',[0.,55.]),
    }
    pressure = pv.read(case/'field_diagnostics/data/pressure_surface_si.vtp')
    wall = pv.read(case/'field_diagnostics/data/wall_wss_si.vtp')
    field = pv.read(case/'frozen_flow/steady_flow_mean_2p0_mmps.vtu')
    ids = pressure['GlobalNodeID_zero_based']
    assert np.allclose(pressure.points,field.points[ids],rtol=0,atol=2e-11)
    assert np.array_equal(pressure['Pressure_Pa'],field['Pressure'][ids])
    faces = wall.faces.reshape(-1,4)[:,1:]
    area = wall.cell_data['Area_m2']
    raw = wall.cell_data['WSS_raw_Pa']
    traction = wall.cell_data['Tangential_viscous_traction_Pa']
    assert np.allclose(np.linalg.norm(traction,axis=1),raw,rtol=1e-13,atol=1e-13)
    assert np.max(np.abs(np.einsum('ij,ij->i',traction,wall['Outward_normal'])))<1e-10
    weights = np.bincount(faces.ravel(),weights=np.repeat(area,3),minlength=wall.n_points)
    total = np.bincount(faces.ravel(),weights=np.repeat(raw*area,3),minlength=wall.n_points)
    reconstructed = total/weights
    assert np.allclose(reconstructed,wall['WSS_display_Pa'],rtol=1e-13,atol=1e-13)
    cameras,annotations = {},{}
    for kind,surface in [('pressure',pressure),('wss',wall)]:
        record = media['views'][kind]
        footer,unit,limits = expected[kind]
        assert record['footer']==footer and record['colorbar_title']==unit
        assert record['colorbar_range_Pa']==limits and not record['colorbar_extra_text']
        assert record['magnification_vs_previous_pressure_wss']>1.3
        trace = json.loads((output/f'{kind}_camera.json').read_text())
        positions = np.array([r['camera_position_um'] for r in trace])
        centers = np.array([r['center_um'] for r in trace])
        up = np.array([r['camera_up_unit'] for r in trace])
        scales = np.array([r['parallel_scale_um'] for r in trace])
        assert len(trace)==432 and np.all(up==[0,0,1])
        assert np.all(centers==centers[:1]) and np.ptp(scales)==0.
        assert np.ptp(positions[:,2])<1e-12
        assert np.ptp(np.linalg.norm((positions-centers)[:,:2],axis=1))<1e-10
        steps = np.diff(np.r_[[r['angle_deg'] for r in trace],trace[0]['angle_deg']+360])
        assert np.allclose(steps,360/432)
        # Reconstruct projection from the actual camera, independent of the
        # renderer's recorded bounds, and check every surface point every frame.
        maximum_bound = 0.
        for cam in trace:
            direction = np.array(cam['camera_position_um'])-cam['center_um']
            direction /= np.linalg.norm(direction)
            right = np.cross([0.,0.,1.],direction);right /= np.linalg.norm(right)
            vertical = np.cross(direction,right)
            q = np.asarray(surface.points,dtype=float)*1e6-cam['center_um']
            assert abs(right[2])<1e-14 and vertical[2]>0.
            x = q@right/(cam['parallel_scale_um']*(1920*.90/(1080*.86)))
            y = q@vertical/cam['parallel_scale_um']
            bound = float(max(np.abs(x).max(),np.abs(y).max()))
            assert bound<=.975001
            maximum_bound = max(maximum_bound,bound)
        cameras[kind] = dict(z_axis_vertical=True,fixed_axis=[0,0,1],
            fixed_center_height_radius_scale=True,all_surface_points_inside_view=True,
            max_normalized_projected_bound=maximum_bound,
            magnification_vs_previous=record['magnification_vs_previous_pressure_wss'])
        annotations[kind] = validate_annotations(case,output,kind,trace,record['source_surface'])
    assert (output/'pressure_camera.json').read_bytes()==(output/'wss_camera.json').read_bytes()
    for video in media['videos']:
        assert video['decoded_frames']==video['distinct_frames']==432
        assert video['size']==[1920,1080] and video['fps']==24 and video['duration_s']==18.
        assert video['no_extra_colorbar_text_all_frames']
        assert sha(output/video['file'])==video['sha256']
    result = dict(all_pass=True,source_files_unchanged=True,
        pressure_matches_frozen_FEM_pressure_exactly=True,
        wss_nodal_area_average_max_error_Pa=float(np.abs(reconstructed-wall['WSS_display_Pa']).max()),
        wss_raw_magnitudes_and_tangential_tractions_consistent=True,
        camera_checks=cameras,annotation_checks=annotations,
        requested_footer_strings_exact=True,colorbars_only_titles_and_numeric_ticks=True,
        identical_camera_and_scale_for_pressure_and_wss=True,
        videos_hash_verified=True,renderer_and_shared_code_hash_verified=True)
    (output/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,default=HERE/'input_data')
    parser.add_argument('--output',type=Path,default=HERE/'results/surface_fields')
    args = parser.parse_args()
    validate(args.case.resolve(),args.output.resolve())
