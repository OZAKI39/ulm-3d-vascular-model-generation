from pathlib import Path
import numpy as np
import pytest
import pyvista as pv
from particle_3d.particle8_replay import REPO,read,digest,snapshot,interpolate
from particle_3d.particle8_full3d_data import (CASES,verify_preservation,discover,trajectory_tail,schedule,camera_angles,replay_frame,safe_output)
from particle_3d.particle8_full3d_visuals import display_position,shape_mesh


def test_all_upstream_files_still_identical(root):assert verify_preservation(root)>500


def test_discovered_sources_are_authoritative(root):
    assert discover()==read(root/'data/discovered_inputs.json')['scenes']
    assert len(discover())==4


def test_complete_frozen_wall_and_boundary_geometry(root):
    data=np.load(root/'data/full_frozen_geometry.npz');m=read(root/'data/full_geometry_manifest.json')
    assert m['full_wall_triangle_count']==45221 and not m['decimation']
    for role,spec in m['surfaces'].items():
        path=REPO/'formal_3D_flow_solver/FEM_SimVascular'/spec['source']['path']
        assert digest(path)==spec['source']['sha256'];s=pv.read(path)
        np.testing.assert_array_equal(data[role+'_points_m'],s.points)
        np.testing.assert_array_equal(data[role+'_faces'],s.faces.reshape(-1,4)[:,1:])


def test_recent_tail_includes_only_saved_knots_and_clipped_endpoints(scenes):
    r=scenes['real_single_mb_entry']['records'][0];t=r['trajectory'][-1]['time_s'];window=.00035
    tail=trajectory_tail(r,t,window);assert tail[0]['time_s']==t-window and tail[-1]['time_s']==t
    saved={s['time_s'] for s in r['trajectory']}
    assert all(s['time_s'] in saved for s in tail[1:-1])
    for s in tail:assert {k:s[k] for k in ['position_m','q']}==interpolate(r,s['time_s'])
    full=trajectory_tail(r,t,None);assert full==r['trajectory']


def test_pending_and_deleted_never_have_visible_tail(scenes):
    r=scenes['real_mixed_inlet_smoke']['records'][-1];assert not trajectory_tail(r,42.986,.001)
    r=scenes['synthetic_open_section']['records'][0]
    assert not trajectory_tail(r,r['scheduled_time_s']-1e-8,.001)
    assert not trajectory_tail(r,r['delete_time_s'],.001)


@pytest.mark.parametrize('window',[0,-1,float('nan'),float('inf')])
def test_bad_tail_windows_rejected(scenes,window):
    r=scenes['real_single_mb_entry']['records'][0]
    with pytest.raises(ValueError):trajectory_tail(r,r['trajectory'][-1]['time_s'],window)


def test_static_rbc_does_not_gain_motion_from_tail(scenes):
    scene=scenes['real_mixed_inlet_smoke'];r=scene['records'][0]
    for t in [1.,20.,42.986]:
        tail=trajectory_tail(r,t,1.)
        assert all(s['position_m']==r['trajectory'][0]['position_m'] for s in tail)
        frame=replay_frame(scene,t,1.)
        assert not set(frame['pending_ids'])&{p['particle_id'] for p in frame['particles']}
        assert frame['pending_drawn_ids']==[]


@pytest.mark.parametrize('case',list(CASES))
def test_monotone_camera_and_physical_clock(scenes,case):
    times=schedule(case,scenes[CASES[case]['scene']]);t=np.array([x[0] for x in times]);assert np.all(np.diff(t)>=0)
    angles=np.array([camera_angles(case,i,len(times))['azimuth_deg'] for i in range(len(times))])
    assert angles[0]==35
    if CASES[case]['camera']=='fixed':assert np.all(angles==35)
    else:assert angles[-1]==115 and np.all(np.diff(angles)>0) and np.ptp(np.diff(angles))<1e-12


def test_display_transform_is_separate_from_physical_coordinates():
    x=np.array([[1e-5,2e-5,3e-5]]);original=x.copy();origin=np.array([0,0,1e-5])
    actual=display_position(x,origin);stretched=display_position(x,origin,True)
    np.testing.assert_array_equal(x,original);np.testing.assert_allclose(actual,[ [10,20,20] ])
    np.testing.assert_allclose(stretched,actual*[1,1,20000])


def test_original_shapes_and_quaternions_are_not_mutated(scenes):
    from particle_3d.particle8_replay import canonical_hash
    from particle_3d.rbc_orientation import rotation_matrix
    for name,t in [('real_single_mb_entry',42.986),('real_mixed_inlet_smoke',42.),('synthetic_open_section',.11785)]:
        scene=scenes[name];before=canonical_hash(scene);records={r['particle_id']:r for r in scene['records']}
        for p in snapshot(scene,t)['active']:
            record=records[p['particle_id']];shape=record['admitted_shape']
            mesh=shape_mesh(record,p,np.zeros(3))
            local=np.asarray(mesh.points,dtype=float)-np.array(p['position_m'])*1e6
            # A tessellated, rotated capsule's bounding-box center is not its
            # physical center. Verify the analytic surface about the saved
            # center and axis/orientation instead of the mesh bounding box.
            if shape['mode']=='SPHERE_MB':
                norm=np.linalg.norm(local,axis=1)/(shape['radius_m']*1e6)
            elif shape['mode']=='FREE_OBLATE':
                body=local@rotation_matrix(p['q'])
                norm=np.linalg.norm(body/(np.array(shape['axes_m'])*1e6),axis=1)
            else:
                axis=np.array(shape['capsule_axis']);half=shape['capsule_length']*1e6/2
                axial=np.clip(local@axis,-half,half)
                norm=np.linalg.norm(local-axial[:,None]*axis,axis=1)/(shape['capsule_radius']*1e6)
            np.testing.assert_allclose(norm,1.,atol=1e-5,rtol=0)
        assert canonical_hash(scene)==before


def test_output_guard_prevents_overwriting_old_results():
    for p in [REPO,REPO/'particle_3d/reports',REPO/'particle_3d/reports/particle8/new']:
        with pytest.raises(ValueError):safe_output(p)


def test_real_prebirth_caption_never_claims_synthetic_section(scenes):
    from particle_3d.particle8_full3d_visuals import context_note
    for case in ['A_real_single_mb_orbit','A_real_single_mb_fixed','B_real_mixed_smoke_orbit']:
        scene=scenes[CASES[case]['scene']];state=replay_frame(scene,scene['time_range_s'][0],.00035)
        note=context_note(case,state,scene)
        assert not state['particles'] and '10 nm' not in note and 'straddle' not in note
        assert 'no active' in note or 'No admitted' in note
