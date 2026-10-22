"""Geometric positive and negative controls for the compact C stage."""
from pathlib import Path
import numpy as np
import cadquery as cq
import networkx as nx
from .support_inputs import ROOT,save
from .casting_fixture_geometry import bounds_box,shape_mesh
from .casting_fixture_v2_geometry import world
from .casting_mold_qc import fcl_distance
from .support_removal import swept_path
from .compact_c_geometry import local_cap,validate_parameters
from .all_abs_geometry import conformal_land
from .all_abs_audit import contact_surface_evidence

NAMES=['single_port','master_port','floating_secondary','too_tight_C','closed_ring_invalid','C_hits_tree','C_insertion_blocked','wall_solvent_route_pass','wall_solvent_route_blocked']


def run_case(name,out):
    if name not in NAMES:raise ValueError(name)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    p=dict(port_id='I1',role='MASTER' if name=='master_port' else 'SECONDARY',radius_mm=1.,anti_float_gap_mm=.22,cradle_clearance_mm=.33,wall_distance_from_trim_mm=6.,axis=np.array([0.,1.,0.]),trim_origin=np.zeros(3))
    cfg=dict(retainer_wall_thickness_lines=2,retainer_axial_length_mm=1.2,retainer_cap_front_from_trim_mm=None,retainer_neck_overlap_mm=.2,retainer_neck_ramp_length_mm=1.3,retainer_neck_half_width_mm=.3,retainer_wallward_slide_mm=1.6)
    old=dict(design=[dict(port_id='I1',retainer_guide_clearance_mm=.12)],local_parts={'I1':dict(retainer=bounds_box([-1,1.05,1.22],[1,5.76,1.88]),inner_guide_x=1.88)})
    angle=360 if name=='closed_ring_invalid' else 90;clear=-.02 if name=='too_tight_C' else .22;cfg['port_clearance_override_mm']=clear;error=None
    try:cap,params=local_cap(p,old,cfg,.22,angle)
    except ValueError as exc:error=str(exc);cfg['port_clearance_override_mm']=.22;cap,params=local_cap(p,old,cfg,.22,90)
    core=cq.Solid.makeCylinder(1.,9.,cq.Vector(0,-4.5,0),cq.Vector(0,1,0));cm=shape_mesh(core);capm=shape_mesh(cap);tree=bounds_box([-.3,1.4,1.2],[.3,1.7,1.9]) if name=='C_hits_tree' else bounds_box([-.5,-5,-.5],[.5,-4,.5]);tm=shape_mesh(tree)
    start=np.eye(4);start[2,3]=5.;finish=np.eye(4);obstacles=[cm];block=bounds_box([-2,1,3],[2,6.2,3.25])
    if name=='C_insertion_blocked':obstacles.append(shape_mesh(block))
    path=swept_path(capm,obstacles,[start,finish]);up=cm.copy();up.apply_translation([0,0,1]);upper=fcl_distance(up,capm)[0]<=0
    wall=bounds_box([-2,6,-3],[2,8,4]);beam=bounds_box([-.33,1.7,-3],[.33,6.2,-1.4]);local_land,record=conformal_land(cm,p,.66,1.,2.2);land=world(local_land,p)
    if name=='wall_solvent_route_blocked':land=land.translate((0,0,-.2))
    ev=contact_surface_evidence(cm,shape_mesh(land),record['contact_witnesses'],record['contact_area_mm2'],dict(contact_area_tolerance_mm2=.01,contact_distance_tolerance_mm=1e-5));G=nx.DiGraph();G.add_nodes_from(['bath','wall','blade','contact','port']);G.add_edge('bath','wall')
    if wall.intersect(beam).Volume()>0:G.add_edge('wall','blade')
    if beam.intersect(land).Volume()>0:G.add_edge('blade','contact')
    if ev['passed']:G.add_edge('contact','port')
    # Test the approved +/-1 fit interval against a real end stop only for master.
    stop=bounds_box([-.4,4.72,-3],[.4,6,.5]);ax=[]
    for delta in [-1.,1.]:
        probe=cm.copy();probe.apply_translation([0,delta,0]);ds=[fcl_distance(probe,capm)[0],fcl_distance(probe,shape_mesh(wall))[0]]
        if p['role']=='MASTER':ds.append(fcl_distance(probe,shape_mesh(stop))[0])
        ax.append(dict(delta_mm=delta,clear=min(ds)>0))
    for k,part in dict(core=core,C_retainer=cap,central_tree=tree,wall=wall,blade=beam,contact_land=land,insertion_blocker=block,master_stop=stop).items():shape_mesh(part).export(out/(k+'.stl'))
    if name=='closed_ring_invalid':shape_mesh(cq.Workplane('XZ').circle(1.66).circle(1.22).extrude(1.2).val()).export(out/'rejected_closed_ring.stl')
    result=dict(case=name,parameter_rejected=error is not None,rejection_reason=error,parameters=params,positive_core_gap_mm=fcl_distance(cm,capm)[0],central_tree_contact=fcl_distance(tm,capm)[0]<=0,insertion=path,no_snap_required=True,port_up_1mm_blocked=upper,role=p['role'],axial_probes=ax,secondary_float_pass=all(r['clear'] for r in ax),master_positive_stop=(not ax[1]['clear']),wall_route_pass=nx.has_path(G,'bath','port'),contact_evidence=ev,graph_nodes=list(G),graph_edges=list(G.edges),no_retainer_route=True)
    save(out/'input.json',dict(case=name,wrap_angle_deg=angle,upper_clearance_mm=clear,radius_mm=1.,scope='Synthetic dimensions only; not anatomy'));save(out/'result.json',result);return result


def generate(out=None):
    out=Path(out or ROOT/'tests/data/all_abs_casting/c_retainer');return [run_case(n,out/n) for n in NAMES]
