"""Two-top-port alternative: one smooth U-shaped centerline, no horizontal stub.

Keep the established inlet and OUTLET1 axes and the low microscopy segment.
Replace the tee/OUTLET2 by an upward circular elbow at the same OUTLET1 axis.
The frame, raised dovetails, native slicing and Rhino export are shared.
"""
from __future__ import annotations

import cadquery as cq
import numpy as np
import trimesh

from . import y_channel_casting as shared


def design(cfg):
    r=cfg['diameter_mm']/2;R=r*cfg['bend_radius_to_lumen_radius']
    z=cfg['base_mm']+cfg['pdms_floor_mm']+r
    start=-cfg['parent_length_mm'];outlet=cfg['downstream_length_mm']
    left=start-cfg['terminal_straight_mm'];inlet=left-R;right=outlet-R
    if R<2*r or right<=left:
        raise ValueError('TWO_PORT_ELBOW_RADIUS_OR_STRAIGHT_LENGTH_INVALID')
    inner=np.array([[inlet-r-cfg['side_clearance_mm'],-r-cfg['side_clearance_mm']],
                    [outlet+r+cfg['side_clearance_mm'],r+cfg['side_clearance_mm']]])
    inner=shared.enforce_pdms_width(inner,cfg)
    outer=inner+[[-cfg['wall_mm']]*2,[cfg['wall_mm']]*2]
    fill=max(z+R+cfg['pdms_above_elbow_mm'],cfg['base_mm']+shared.minimum_cad(cfg,'height'))
    rim=fill+cfg['freeboard_mm'];collar_z=shared.retainer_arm_bottom(cfg,rim)+cfg['retainer_thickness_mm']
    return dict(radius=r,bend_radius=R,center_z=z,start_x=start,junction_x=outlet,
        inlet_x=inlet,bend_x=left,outlet_bend_x=right,straight_centerline_length_mm=right-left,
        inner=inner,outer=outer,fill=fill,rim=rim,collar_z=collar_z,top=collar_z+cfg['port_above_retainer_mm'],
        gate_back=None,horizontal_end=None,chamber_mode=cfg['chamber_mode'],
        channel_topology='two_vertical_ports',horizontal_outlet_present=False)


def build_core(cfg,d):
    r,R,z=d['radius'],d['bend_radius'],d['center_z'];top=d['top']
    v=lambda x,h:cq.Vector(x,0,h)
    inlet_top=v(d['inlet_x'],top);inlet_knee=v(d['inlet_x'],z+R)
    inlet_mid=v(d['bend_x']-R/np.sqrt(2),z+R*(1-1/np.sqrt(2)))
    inlet_low=v(d['bend_x'],z);outlet_low=v(d['outlet_bend_x'],z)
    outlet_mid=v(d['outlet_bend_x']+R/np.sqrt(2),z+R*(1-1/np.sqrt(2)))
    outlet_knee=v(d['junction_x'],z+R);outlet_top=v(d['junction_x'],top)
    wire=cq.Wire.assembleEdges([cq.Edge.makeLine(inlet_top,inlet_knee),
        cq.Edge.makeThreePointArc(inlet_knee,inlet_mid,inlet_low),
        cq.Edge.makeLine(inlet_low,outlet_low),
        cq.Edge.makeThreePointArc(outlet_low,outlet_mid,outlet_knee),
        cq.Edge.makeLine(outlet_knee,outlet_top)])
    channel=cq.Solid.sweep(cq.Wire.makeCircle(r,inlet_top,cq.Vector(0,0,-1)),[],wire).clean()
    if not channel.isValid() or len(channel.Solids())!=1:
        raise ValueError('INVALID_TWO_TOP_PORT_SWEEP')
    core=channel
    for x in [d['inlet_x'],d['junction_x']]:
        core=core.fuse(cq.Solid.makeCylinder(r+cfg['collar_radial_width_mm'],cfg['collar_thickness_mm'],
            v(x,d['collar_z']),cq.Vector(0,0,1))).clean()
    return core,channel


def section_samples(d):
    r,R,z=d['radius'],d['bend_radius'],d['center_z']
    rows=[dict(region='straight',center_mm=[x,0,z],normal=[1,0,0])
          for x in np.linspace(d['bend_x']+.5,d['outlet_bend_x']-.5,7)]
    for side,sign,base in [('INLET',-1,d['bend_x']),('OUTLET1',1,d['outlet_bend_x'])]:
        for angle in np.linspace(.12,np.pi/2-.12,5):
            rows.append(dict(region=side+'_elbow',center_mm=[base+sign*R*np.sin(angle),0,z+R*(1-np.cos(angle))],
                normal=[sign*np.cos(angle),0,np.sin(angle)]))
        x=d['inlet_x'] if side=='INLET' else d['junction_x']
        for h in [z+R+.5,d['fill']-.5,d['top']-.4]:
            rows.append(dict(region=side+'_riser',center_mm=[x,0,h],normal=[0,0,1]))
    return rows


def centerline_samples(d):
    """Dense analytical elbow/straight samples used against actual extrusion."""
    R,z=d['bend_radius'],d['center_z']
    rows=[[x,0,z] for x in np.linspace(d['bend_x'],d['outlet_bend_x'],80)]
    for sign,base in [(-1,d['bend_x']),(1,d['outlet_bend_x'])]:
        rows.extend([[base+sign*R*np.sin(a),0,z+R*(1-np.cos(a))] for a in np.linspace(0,np.pi/2,40)])
    return rows


def channel_qc(d,channel,ports):
    """Independent topology/volume checks on the tessellated lumen surface."""
    top=np.max(abs(channel.triangles[:,:,2]-d['top']),axis=1)<1e-5
    wall=channel.copy();wall.update_faces(~top);wall.remove_unreferenced_vertices()
    edge_count=np.bincount(wall.edges_unique_inverse)
    boundary=wall.edges_unique[edge_count==1]
    loops=trimesh.graph.connected_components(boundary,min_len=1)
    degrees=np.bincount(boundary.ravel()) if len(boundary) else np.array([])
    limit=d['junction_x']+d['radius']
    expected_length=d['straight_centerline_length_mm']+np.pi*d['bend_radius']+2*(d['top']-d['center_z']-d['bend_radius'])
    expected_volume=np.pi*d['radius']**2*expected_length
    relative_volume_error=abs(float(channel.volume)-expected_volume)/expected_volume
    passed=(len(ports)==2 and {p['name'] for p in ports}=={'INLET','OUTLET1'}
            and all(p['outward_axis']==[0,0,1] for p in ports)
            and len(loops)==2 and np.all(degrees[degrees>0]==2)
            and channel.bounds[1,0]<=limit+1e-5 and relative_volume_error<.003)
    return dict(passed=bool(passed),top_opening_count=len(loops),horizontal_outlet_present=False,
        max_channel_x_mm=float(channel.bounds[1,0]),expected_max_channel_x_mm=limit,
        straight_centerline_length_mm=d['straight_centerline_length_mm'],
        analytic_centerline_length_mm=expected_length,analytic_channel_volume_mm3=expected_volume,
        mesh_channel_volume_mm3=float(channel.volume),volume_relative_error=relative_volume_error,
        construction='Single circular sweep with two tangent 90-degree elbows; no tee or blind stub',
        manual_horizontal_cut_required=False)
