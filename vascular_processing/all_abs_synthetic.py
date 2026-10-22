"""Small geometric counterexamples for solvent access and peripheral trimming.

Rectangular lumen surrogates make finite-area contact exact, without invoking
any vascular fitting or production source. Fixtures and core remain separate.
"""
import json
from pathlib import Path
import numpy as np
import networkx as nx
from .support_inputs import ROOT,save
from .casting_fixture_geometry import bounds_box,shape_mesh
from .casting_fixture_v2_audit import cad_manifold
from .casting_mold_qc import fcl_distance
from .all_abs_audit import graph_status,network_access

NAMES=['single_port_access','two_port_access','branch_network_access','fixture_only_on_extension','fixture_touches_tree_fail','blocked_port_access','trim_removes_access_void','trim_hits_branch_fail','vapor_smoothing_scope']


def vapor_scope(requested):
    allowed={'INNER_BOTTOM','INNER_SIDE_WALLS','OUTER_MOLD'}
    return dict(passed=set(requested)<=allowed,core_forbidden=True,status='VAPOR_SMOOTHING_NOT_DIMENSIONALLY_CALIBRATED')


def run_case(name,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    core=bounds_box([-3,-.5,2.5],[3,.5,3.5]);central=bounds_box([-.8,-.5,2.5],[.8,.5,3.5]);wall=bounds_box([-6,-2,0],[-5,2,5]);beam=bounds_box([-5.2,-.2,0],[-2,.2,2.0 if name=='blocked_port_access' else 2.5]);beam2=bounds_box([2,-.2,0],[5.2,.2,2.5]);rightwall=bounds_box([5,-2,0],[6,2,5])
    if name=='branch_network_access':
        branch=bounds_box([-.25,0,2.75],[.25,1.5,3.25]);core=core.fuse(branch);central=central.fuse(branch)
    if name=='fixture_touches_tree_fail':beam=bounds_box([-.5,-.2,0],[.5,.2,2.5])
    cm=cad_manifold(core);bm=cad_manifold(beam);eps=1e-4;intersection=float((bm^cm).volume());gap=fcl_distance(shape_mesh(core),shape_mesh(beam))[0];contact=float((cad_manifold(beam.translate((0,0,eps)))^cm).volume())/eps
    G=nx.DiGraph();G.add_node('external_bath');ports=['port_I1'];G.add_node('port_I1');G.add_edge('external_bath','wall_I1')
    root=float(beam.intersect(wall).Volume())
    if root>1e-7:G.add_edge('wall_I1','fixture_I1')
    if gap<1e-7 and contact>.01 and intersection<1e-7:G.add_edge('fixture_I1','port_I1')
    if name in ['two_port_access','branch_network_access']:
        ports.append('port_O1');G.add_node('port_O1');G.add_edge('external_bath','wall_O1')
        rightroot=float(beam2.intersect(rightwall).Volume());rightgap=fcl_distance(shape_mesh(core),shape_mesh(beam2))[0]
        rightarea=float((cad_manifold(beam2.translate((0,0,eps)))^cm).volume())/eps
        if rightroot>1e-7:G.add_edge('wall_O1','fixture_O1')
        if rightgap<1e-7 and rightarea>.01 and (cad_manifold(beam2)^cm).volume()<1e-7:G.add_edge('fixture_O1','port_O1')
        shape_mesh(beam2).export(out/'fixture_O1.stl');shape_mesh(rightwall).export(out/'wall_O1.stl')
    access=graph_status(G,ports);treegap=fcl_distance(shape_mesh(central),shape_mesh(beam))[0];cut=bounds_box([-6,-2,-1],[.2 if name=='trim_hits_branch_fail' else -1.8,2,6]);cutm=cad_manifold(cut);remaining=float((bm-cutm).volume());hits=float((cutm^cad_manifold(central)).volume())
    data=np.array([[1,1,-3,0,3,.5,-1],[2,1,0,0,3,.5,1],[3,1,3,0,3,.5,2],[4,1,0,1,3,.3,2]])
    np.savetxt(out/'network.swc',data,fmt=['%d','%d','%.6f','%.6f','%.6f','%.6f','%d']);bad=np.r_[data,[[99,1,12,0,3,.2,-1]]];np.savetxt(out/'disconnected.swc',bad,fmt=['%d','%d','%.6f','%.6f','%.6f','%.6f','%d'])
    anchors=[dict(port_id='I1',swc_id=1),dict(port_id='O1',swc_id=3)];reachable=[x.replace('port_','') for x in access['reachable_ports']];network=network_access(out/'network.swc',anchors,reachable);disconnected=network_access(out/'disconnected.swc',anchors,reachable)
    for key,part in dict(core=core,central_tree=central,fixture=beam,wall=wall,trim_region=cut).items():shape_mesh(part).export(out/(key+'.stl'))
    result=dict(case=name,nominal_geometry='Rectangular-lumen synthetic surrogate; NOT production anatomy',solvent_status=access['status'],reachable_ports=access['reachable_ports'],wall_attachment_volume_mm3=root,core_gap_mm=gap,contact_area_mm2=contact,interference_mm3=intersection,central_tree_contact=treegap<=1e-7,auxiliary_after_trim_mm3=remaining,central_cut_volume_mm3=hits,branch_network=network,disconnected_network_rejected=not disconnected['all_vascular_branches_reachable'],allowed_smoothing=vapor_scope(['INNER_BOTTOM','INNER_SIDE_WALLS']),core_smoothing_rejected=not vapor_scope(['VASCULAR_CORE'])['passed'])
    save(out/'input.json',dict(case=name,blocked_gap_mm=.5 if name=='blocked_port_access' else 0,expected_access='SOLVENT_ACCESS_BLOCKED' if name in ['blocked_port_access','fixture_touches_tree_fail'] else 'SOLVENT_ACCESS_GEOMETRIC_PASS'));save(out/'result.json',result);return result


def generate(out=None):
    out=Path(out or ROOT/'tests/data/all_abs_casting');return [run_case(n,out/n) for n in NAMES]
