"""V2 profile and casting-frame inputs; V1 detection remains reusable."""
import copy,json
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import ROOT,discover,resolve_profile,transform,swc_branches,port_branches,save,numeric,digest

def calibration_signature(profile):
    # Bridge performance changes with speed, flow and temperature as well as
    # nozzle size. An unrelated ABS coupon must not silently unlock a span.
    return {k:profile[k] for k in ['printer','process_name','nozzle_diameter','layer_height','line_width','material','bridge_settings','speeds','temperatures']}

def calibrated_span(calibration,profile):
    if calibration.get('physical_test_passed') is not True or calibration.get('profile')!=calibration_signature(profile):return None
    span=float(calibration.get('maximum_bridge_span_mm',0.))
    return span if np.isfinite(span) and span>0 else None

def load(out,source=None):
    inputs=discover(source=source);profile,discovery=resolve_profile(inputs,out/'profiles')
    m=discovery['flattened']['machine'];p=discovery['flattened']['process'];f=discovery['flattened']['filament']
    xy=np.array([[float(x) for x in v.split('x')] for v in m['printable_area']])
    profile['build_bounds_mm']=[[float(xy[:,0].min()),float(xy[:,1].min()),0.],[float(xy[:,0].max()),float(xy[:,1].max()),numeric(m['printable_height'])]]
    profile['bridge_settings']={k:v for k,v in p.items() if 'bridge' in k}
    profile['speeds']={k:p.get(k) for k in ['outer_wall_speed','inner_wall_speed','bridge_speed']}
    profile['temperatures']={k:f.get(k) for k in ['nozzle_temperature','textured_plate_temp']}
    calibration=ROOT/'config/fdm_bridge_calibration.json';profile['bridge_limit_mm']=2.;profile['bridge_calibration_status']='BRIDGE_LIMIT_UNCALIBRATED'
    if calibration.is_file():
        c=json.loads(calibration.read_text())
        limit=calibrated_span(c,profile)
        if limit is not None:
            profile['bridge_limit_mm']=limit;profile['bridge_calibration_status']='PHYSICAL_COUPON_CALIBRATED'
    profile['gap_status']='PHYSICAL_GAP_UNVALIDATED';profile['roof_gap_mm']=.2
    mold=trimesh.load_mesh(inputs['mold']);vascular=trimesh.load_mesh(inputs['core']);box=trimesh.load_mesh(inputs['box'])
    casting=np.array(json.loads(inputs['casting_transform'].read_text())['transform_4x4'])
    branches,topology=swc_branches(inputs['smoothed_swc'],casting,inputs['branch_mapping'],.25)
    inv=np.linalg.inv(inputs['print_transform'])
    for b in port_branches(inputs,.25):
        b['xyz']=transform(b['xyz'],inv);b['tangent']=np.asarray(b['tangent'])@inv[:3,:3].T;branches.append(b)
    save(out/'profile_runtime.json',profile);save(out/'source_discovery.json',inputs)
    save(out/'v2_policy_runtime.json',dict(architecture='FDM_NATIVE_VASCULAR_SUPPORT_V2',scale_factor=1.,casting_orientation_frozen=True,printing_orientation_fixed=False,
        root_domains=['BUILD_PLATE','BOX_BOTTOM','BOX_SIDE_WALL'],support_count_limit=None,bifurcation_policy='DAMAGE_REMOVAL_PENALTY_NOT_HARD_KEEPOUT',
        primary_gate='COMPLETE_LAYER_POLYGON_SUPPORTABILITY',legacy_5mm_role='COMPARISON_ONLY',legacy_custom_support_config_active=False,
        bridge_calibration_signature=calibration_signature(profile)))
    return inputs,profile,discovery,mold,vascular,box,branches,topology

def oriented(branches,topology,matrix):
    result=copy.deepcopy(branches);top=copy.deepcopy(topology)
    for b in result:
        b['xyz']=transform(b['xyz'],matrix);b['tangent']=np.asarray(b['tangent'])@np.asarray(matrix)[:3,:3].T
    for k in ['bifurcations','endpoint_xyz']:
        if len(top[k]):top[k]=transform(top[k],matrix)
    return result,top

def public(value):
    if isinstance(value,dict):return {k:public(v) for k,v in value.items() if not k.startswith('_')}
    if isinstance(value,list):return [public(v) for v in value]
    return value
