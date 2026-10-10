"""Display-only original native CFD streamlines, terminated at new solid."""
import json
import numpy as np
import pyvista as pv
from .config import dump
from .y_bifurcation_sources import OUT,RESULTS,YSources
from .y_mural_clot_geometry import locator,hit,inside


def prepare_background():
    s=YSources();solid=pv.read(RESULTS/'paraview/clot_solid.vtp');loc=locator(solid)
    flow=s.flow.copy();flow['speed_m_s']=np.linalg.norm(s.u,axis=1)
    flow.field_data['domain_identity']=['UNOBSTRUCTED_STEADY_BACKGROUND_FLOW']
    flow.field_data['CLOT_PRESENT_CFD_SOLVED']=['NOT_RUN']
    flow.field_data['warning']=['Original unobstructed background CFD; not recomputed with clot']
    flow.save(RESULTS/'paraview/background_flow.vtu')
    arrow_positions=np.array([[-.010,0,0],[.005,.005*np.tan(np.pi/6),0],[.005,-.005*np.tan(np.pi/6),0]])
    arrow_velocity,_=s.sample(arrow_positions)
    np.savez_compressed(RESULTS/'FIGURE_DATA/flow_direction_arrows.npz',position_xyz_m=arrow_positions,velocity_m_s=arrow_velocity)
    # Fixed modest inlet distribution, plus near +z wall seeds to expose clipping.
    seed=[]
    for radius,count in [(.00025,8),(.00055,12),(.00085,16),(.00115,20)]:
        for theta in np.linspace(0,2*np.pi,count,endpoint=False):seed.append([-.01375,radius*np.sin(theta),radius*np.cos(theta)])
    seed=np.array(seed);s.sample(seed)  # strict membership before tracing
    raw=flow.streamlines_from_source(pv.PolyData(seed),vectors='Velocity',integration_direction='forward',
         initial_step_length=5e-5,min_step_length=1e-5,max_step_length=1e-4,step_unit='l',max_length=.060,
         max_steps=20000,terminal_speed=1e-8,compute_vorticity=False)
    raw.save(RESULTS/'FIGURE_DATA/original_native_streamlines.vtp')
    points=[];lines=[];vel=[];truncated=0;pos=0;count=0
    r=raw.lines
    while pos<len(r):
        n=int(r[pos]);ids=r[pos+1:pos+1+n];pos+=n+1
        if len(ids)<2:continue
        path=raw.points[ids];kept=[path[0]]
        for a,b in zip(path[:-1],path[1:]):
            ok,x,_=hit(loc,a,b)
            if ok:
                # Stop immediately before contact; never continue/re-route behind clot.
                d=x-a;L=np.linalg.norm(d)
                if L>1e-9:kept.append(x-d/L*1e-9)
                truncated+=1;break
            kept.append(b)
        if len(kept)<2:continue
        start=len(points);points.extend(kept);lines.extend([len(kept),*range(start,start+len(kept))]);count+=1
    points=np.asarray(points);u,p=s.sample(points)
    visible=pv.PolyData();visible.points=points;visible.lines=np.array(lines)
    visible['speed_m_s']=np.linalg.norm(u,axis=1);visible['velocity_m_s']=u;visible['pressure_Pa']=p
    visible.save(RESULTS/'paraview/background_streamlines_masked.vtp')
    # Endpoint and segment checks against SAME physical clot, no invisible bridge.
    mids=[];pos=0
    while pos<len(visible.lines):
        n=visible.lines[pos];q=points[visible.lines[pos+1:pos+1+n]];pos+=n+1;mids.extend((q[:-1]+q[1:])/2)
    mids=np.asarray(mids)
    assert not inside(solid,mids).any()
    np.savez_compressed(RESULTS/'FIGURE_DATA/flow_display_arrays.npz',position_xyz_m=points,lines=visible.lines,speed_m_s=visible['speed_m_s'],
                        velocity_xyz_m_s=u,pressure_Pa=p,seeds_xyz_m=seed,source_flow_id=s.flow_id)
    report=dict(status='PASS',source_flow_id=s.flow_id,streamline_count=count,seed_count=len(seed),
                streamlines_terminated_at_clot=truncated,visible_segment_midpoints_checked=len(mids),visible_segments_inside_clot=0,
                label='UNOBSTRUCTED BACKGROUND CFD — NOT CLOTTED FLOW FIELD',
                operation='Native VTK streamline integration; first exact intersection with fitted clot terminates displayed line. No detour is constructed.',
                original_field_unchanged=True,velocity_units='m/s; display mm/s',pressure_units='Pa',
                source_velocity_range_m_s=[float(np.linalg.norm(s.u,axis=1).min()),float(np.linalg.norm(s.u,axis=1).max())],
                source_pressure_range_Pa=[float(s.p.min()),float(s.p.max())],CLOT_PRESENT_CFD_SOLVED='NOT_RUN',CLOT_SURFACE_TRACTION_VALIDATED='NOT_ESTABLISHED')
    dump(OUT/'BACKGROUND_FLOW_OVERLAY_CHECK.json',report);print('BACKGROUND_READY',count,truncated,flush=True)
    return report
