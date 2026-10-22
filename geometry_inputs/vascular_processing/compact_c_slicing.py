"""Flat axial printing and actual C lip / guide / aperture G-code checks."""
import json
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from scipy.spatial import cKDTree
from .support_inputs import save,table,numeric
from .casting_fixture_geometry import shape_mesh
from .all_abs_slicing import slice_abs
from .bambu_casting_mold import parse_paths,sample_segments


def pack(g):
    R=Rotation.from_euler('x',90,degrees=True).as_matrix();parts=[];records=[];cursor=0.
    for rec in g['compact_parameters']:
        k=rec['port_id'];m=shape_mesh(g['compact_local'][k]);rot=np.eye(4);rot[:3,:3]=R;m.apply_transform(rot);shift=np.array([cursor,0,0])-m.bounds[0];m.apply_translation(shift);cursor=m.bounds[1,0]+5;parts.append(m);records.append(dict(port_id=k,rotation=R,translation=shift,print_axis='LOCAL_EXTENSION_AXIS_Y_TO_BUILD_Z',bed_face='FLAT_AXIAL_FACE_OF_OPEN_C',no_support_required_by_design=True))
    return trimesh.util.concatenate(parts),records


def preservation(record,g,packing,out):
    if not record.get('slice_passed'):return dict(passed=False,status='RETAINER_FEATURE_DROPPED',reason='SLICE_FAILED')
    out=Path(out);parsed=parse_paths(record['archive']);segments=parsed['segments'][parsed['kind']==0];points=sample_segments(segments,.035);tree=cKDTree(points);M=np.asarray(json.loads((out/'print_transform.json').read_text())['transform_4x4']);settings=record['actual_profile'];layer=numeric(settings['layer_height']);first=numeric(settings['initial_layer_print_height']);width=numeric(settings['line_width']);rows=[];slots=[]
    for rec,pk in zip(g['compact_parameters'],packing):
        R=np.asarray(pk['rotation']);shift=np.asarray(pk['translation']);front=rec['cap_front_y_mm'];last=rec['cap_front_y_mm']+rec['total_axial_extent_mm'];levels=np.arange(first,last-front+.001,layer);yc=front+levels-layer/2;yp=front+levels;Ri=rec['upper_clearance_mm']+next(d['radius_mm'] for d in g['compact_ports'] if d['port_id']==rec['port_id']);t=rec['retainer_wall_thickness_mm'];features={}
        arc=[];theta=np.linspace((180-rec['c_retainer_wrap_angle_deg'])/2+5,(180+rec['c_retainer_wrap_angle_deg'])/2-5,19);theta=np.radians(theta)
        for y,cy in zip(yp,yc):
            if front+.03<cy<rec['cap_back_y_mm']-.03:
                arc.extend(np.c_[(Ri+t/2)*np.cos(theta),np.full(len(theta),y),(Ri+t/2)*np.sin(theta)])
        features['C_upper_lip_and_two_line_wall']=np.asarray(arc);guide=[];neck=[]
        for y,cy in zip(yp,yc):
            if rec['neck_end_y_mm']+.04<cy<last-.04:
                xs=np.linspace(-rec['guide_half_width_mm']+.15,rec['guide_half_width_mm']-.15,15);guide.extend(np.c_[xs,np.full(len(xs),y),np.full(len(xs),(rec['guide_bottom_z_mm']+rec['guide_top_z_mm'])/2)])
            if rec['neck_start_y_mm']+.04<cy<rec['neck_end_y_mm']-.04:
                f=(cy-rec['neck_start_y_mm'])/(rec['neck_end_y_mm']-rec['neck_start_y_mm']);z0=(1-f)*(Ri+.05)+f*rec['guide_bottom_z_mm'];z1=(1-f)*(Ri+t-.05)+f*rec['guide_top_z_mm'];neck.append([0,y,(z0+z1)/2])
        features['short_guide']=np.asarray(guide);features['self_supported_neck']=np.asarray(neck)
        for name,qs in features.items():
            bed=(qs@R.T+shift)@M[:3,:3].T+M[:3,3];ds=tree.query(bed)[0];threshold=width*.8;rows.append(dict(port_id=rec['port_id'],feature=name,query_count=len(qs),missing_query_count=int((ds>threshold).sum()),maximum_distance_mm=float(ds.max()),tolerance_mm=threshold,passed=bool(len(qs) and (ds<=threshold).all())))
        # Reverse both rigid transforms; select actual extrusion in this cap's
        # axial slab. A radial bead envelope checks that the passage is open.
        loc=((points-M[:3,3])@M[:3,:3]-shift)@R
        mask=(loc[:,1]>front+.1)&(loc[:,1]<rec['neck_start_y_mm']-.03)&(np.abs(loc[:,0])<Ri+t+.5)&(loc[:,2]>-Ri-.5)&(loc[:,2]<Ri+t+.5)
        ps=loc[mask];rad=np.linalg.norm(ps[:,[0,2]],axis=1);physical_radius=Ri-rec['upper_clearance_mm'];clear=float(rad.min()-width/2-physical_radius) if len(ps) else -1.;slots.append(dict(port_id=rec['port_id'],feature='open_axial_float_passage',model_extrusion_samples=len(ps),estimated_bead_envelope_radial_clearance_mm=clear,nominal_upper_clearance_mm=rec['upper_clearance_mm'],finite_resolution_tolerance_mm=.05,passed=bool(len(ps) and clear>=rec['upper_clearance_mm']-.05)))
    passed=all(x['passed'] for x in rows+slots);result=dict(passed=passed,status='C_RETAINER_FEATURES_PRESERVED' if passed else 'RETAINER_FEATURE_DROPPED',rows=rows,axial_float_slots=slots,model_segment_count=len(segments),support_segment_count=int((parsed['kind']==1).sum()),scope='Every occupied cap/guide layer sampled against actual positive model extrusion; radial aperture includes nominal half-line envelope. Physical bead dimensions uncalibrated.')
    save(out/'C_feature_preservation.json',result);table(out/'C_feature_preservation.csv',rows+slots);return result


def slice_candidate(g,engine,profiles,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);packed,pk=pack(g);packed.export(out/'C_retainers_print_layout.stl');save(out/'packing.json',pk);record=slice_abs(packed,engine,profiles,out,True);features=preservation(record,g,pk,out);return record,features,pk
