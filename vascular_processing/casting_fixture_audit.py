"""Triangle-surface, swept assembly and finite-clearance rigid-body BCF audits."""
from __future__ import annotations
import itertools,json
from pathlib import Path
import numpy as np
import trimesh
from scipy.optimize import linprog
from scipy.spatial.transform import Rotation
from .support_inputs import save,table
from .support_removal import swept_path
from .casting_mold_qc import fcl_distance
from .casting_fixture_geometry import bounds_box,mesh_qc,shape_mesh


def body_constraints(ports):
    pivot=np.mean([p['capture_center'] for p in ports],axis=0);A=[];B=[]
    for p in ports:
        n=np.asarray(p['axis']);u=np.cross(n,[0,0,1.])
        for offset in [-.8,.8]:
            point=np.asarray(p['capture_center'])+n*offset-pivot
            for normal in (u,np.array([0,0,1.])):
                row=np.r_[normal,np.cross(point,normal)]
                A.extend([row,-row]);B.extend([p['radial_clearance_mm']]*2)
        if p['role']=='MASTER_DATUM':
            point=np.asarray(p['end'])-pivot;A.append(np.r_[n,np.cross(point,n)]);B.append(p.get('axial_stop_gap_mm',p['radial_clearance_mm']))
    A=np.asarray(A);B=np.asarray(B);limits={}
    for i,name in enumerate(('X','Y','Z','roll','pitch','yaw')):
        vals=[]
        for sign in (-1,1):
            c=np.zeros(6);c[i]=-sign;r=linprog(c,A_ub=A,b_ub=B,bounds=[(None,None)]*6,method='highs')
            vals.append(float(sign*r.x[i]) if r.success else None)
        limits[name]=dict(negative=vals[0],positive=vals[1],unit='mm' if i<3 else 'rad')
    rank=int(np.linalg.matrix_rank(A));return dict(matrix_rank=rank,limits=limits,pivot=pivot.tolist(),constraint_count=len(B),geometric_Z_lock=limits['Z']['positive'] is not None,pitch_lock=limits['pitch']['positive'] is not None and limits['pitch']['negative'] is not None,roll_lock=limits['roll']['positive'] is not None and limits['roll']['negative'] is not None,physical_force_validated=False,method='Linearized six-DOF halfspaces: radial clearances at both wall faces + unilateral master end stop; secondary axial float is NOT used as a locking constraint'),A,B


def perturbation_fit(ports,errors,floating=True):
    """Deterministic independent port error feasibility with one rigid placement."""
    pivot=np.mean([p['capture_center'] for p in ports],axis=0);A=[];B=[]
    for p,e in zip(ports,errors):
        n=np.asarray(p['axis']);u=np.cross(n,[0,0,1.]);q=np.asarray(p['capture_center'])-pivot
        for v,allow in [(u,p['radial_clearance_mm']),(np.array([0,0,1.]),p['radial_clearance_mm']),(n,p['axial_float_mm'] if floating and p['role']!='MASTER_DATUM' else p['radial_clearance_mm'])]:
            row=np.r_[v,np.cross(q,v)];delta=float(np.asarray(e)@v);A.extend([row,-row]);B.extend([allow-delta,allow+delta])
    return bool(linprog(np.zeros(6),A_ub=A,b_ub=B,bounds=[(None,None)]*6,method='highs').success)


def sweep_translation(mesh,obstacles,vector,tolerance=.004):
    a=np.eye(4);a[:3,3]=vector;b=np.eye(4)
    return swept_path(mesh,obstacles,[a,b],tolerance=tolerance)


def audit(layout,geometry,meshes,cfg,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);frame=meshes['frame'];core=layout['core'];central=layout['central'];parts=[frame,*meshes['closures'].values(),*meshes['pins'].values()];ports=layout['ports'];H=layout['high'][2];lo=layout['low'];hi=layout['high']
    print('BCF audit: core top loading',flush=True)
    top=sweep_translation(core,[frame],[0,0,H-core.bounds[0,2]+5],cfg['sweep_tolerance_mm']);top['direction']=[0,0,-1];top['closure_absent']=True
    contacts=[]
    for name,m in [('FRAME',frame),*[(f'CLOSURE_{k}',v) for k,v in meshes['closures'].items()],*[(f'PIN_{k}',v) for k,v in meshes['pins'].items()]]:
        d,a,b=fcl_distance(central,m);full,ca,cb=fcl_distance(core,m)
        contacts.append(dict(part=name,central_clearance_mm=d,core_clearance_mm=full,central_contact=d<1e-6,closest_core_point=ca.tolist(),closest_PLA_point=cb.tolist()))
    cavity=bounds_box(lo+[1e-6,1e-6,1e-6],hi-[1e-6,1e-6,1e-6]);intrusion=[]
    for name,shape in [('FRAME',geometry['frame']),*geometry['closures'].items(),*[(f'PIN_{k}',v) for k,v in geometry['pins'].items()]]:
        vol=float(shape.intersect(cavity).Volume());intrusion.append(dict(part=name,illegal_intrusion_mm3=vol))
    print('BCF audit: closure and pin sweeps',flush=True)
    sequence=[layout['master_port']]+[p['port_id'] for p in ports if p['port_id']!=layout['master_port']];present=[frame,core];closing=[];locking=[]
    for k in sequence:
        m=meshes['closures'][k];result=sweep_translation(m,present,[0,0,H-m.bounds[0,2]+4],cfg['sweep_tolerance_mm']);closing.append(dict(port_id=k,**result));present.append(m)
        p=next(p for p in ports if p['port_id']==k);pin=meshes['pins'][k];result=sweep_translation(pin,present,np.asarray(p['axis'])*6,cfg['sweep_tolerance_mm']);locking.append(dict(port_id=k,**result));present.append(pin)
    constraints,A,B=body_constraints(ports);manager=trimesh.collision.CollisionManager()
    for i,m in enumerate(parts):manager.add_object(str(i),m)
    finite=[];pivot=np.asarray(constraints['pivot'])
    for name,axis,angle,translation in [('Z+',None,0,[0,0,1.]),('Z-',None,0,[0,0,-1.]),('roll+','x',2,[0,0,0]),('roll-','x',-2,[0,0,0]),('pitch+','y',2,[0,0,0]),('pitch-','y',-2,[0,0,0]),('yaw+','z',2,[0,0,0]),('yaw-','z',-2,[0,0,0]),('X+',None,0,[1,0,0]),('X-',None,0,[-1,0,0]),('Y+',None,0,[0,1,0]),('Y-',None,0,[0,-1,0])]:
        m=np.eye(4)
        if axis:m[:3,:3]=Rotation.from_euler(axis,angle,degrees=True).as_matrix()
        m[:3,3]=pivot-m[:3,:3]@pivot+translation;finite.append(dict(perturbation=name,collision=bool(manager.in_collision_single(core,transform=m))))
    retention=[]
    for k,m in meshes['closures'].items():
        moved=m.copy();moved.apply_translation([0,0,.5]);cm=trimesh.collision.CollisionManager();cm.add_object('FRAME',frame);cm.add_object('PIN',meshes['pins'][k]);retention.append(dict(port_id=k,closure_up_0p5mm_blocked=bool(cm.in_collision_single(moved))))
    constraints['finite_pose_checks']=finite;constraints['closure_lock_checks']=retention;constraints['closures_positive_Z_locked']=all(x['closure_up_0p5mm_blocked'] for x in retention)
    # Entire exact cavity is disjoint from CAD PLA solids, which certifies all
    # cavity bottom points; an independent ray grid catches export mistakes.
    xy=np.array(list(itertools.product(np.linspace(lo[0]+.1,hi[0]-.1,19),np.linspace(lo[1]+.1,hi[1]-.1,19))));origins=np.c_[xy,np.full(len(xy),H/2)]
    allmesh=trimesh.util.concatenate(parts);hits,ids,_=allmesh.ray.intersects_location(origins,np.tile([0,0,-1.],(len(origins),1)),multiple_hits=True);blocked=int(np.count_nonzero(hits[:,2]>=0)) if len(hits) else 0
    cavity_audit=dict(bottom_open=blocked==0 and sum(x['illegal_intrusion_mm3'] for x in intrusion)<1e-8,ray_count=len(xy),blocked_rays=blocked,illegal_fixture_intrusion_volume=sum(x['illegal_intrusion_mm3'] for x in intrusion),individual_intrusions=intrusion,proof='Exact CAD/cavity empty intersection + independent downward mesh rays; no floor/grid/crossbar')
    floor=core.bounds[0,2];wallgap=min(r['central_clearance_mm'] for r in contacts);layout_margin=float(min(central.bounds[0,0]-lo[0],central.bounds[0,1]-lo[1],hi[0]-central.bounds[1,0],hi[1]-central.bounds[1,1]))
    seal_loop=geometry['seal_loop'];seal=dict(bottom_loop_continuous=bool(np.allclose(seal_loop[0],seal_loop[-1])),port_seal_regions=[p['port_id'] for p in ports],physical_leak_tested=False,status='PHYSICAL_SEAL_NOT_VALIDATED',leak_test_status='PHYSICAL_LEAK_TEST_NOT_DONE',channel_depth_mm=cfg['seal_channel_depth_mm'],channel_width_mm=cfg['seal_channel_width_mm'],port_geometry='External split annulus + rectangular moat; mating gaps require externally applied sealant',outside_cavity=True)
    result=dict(top_loading_path=top,closure_insertions=closing,pin_insertions=locking,closure_insertion_pass=all(r['passed'] for r in closing+locking),central_vascular_contact_count=sum(r['central_contact'] for r in contacts),nominal_capture_sites=len(ports),nominal_hard_contact_count=sum(r['core_clearance_mm']<1e-6 for r in contacts),contact_policy='Finite-clearance capture at frozen straight extensions only',contacts=contacts,central_min_clearance_mm=wallgap,central_plane_wall_margin_mm=layout_margin,core_to_dish_clearance_mm=float(floor),anti_float=constraints,cavity=cavity_audit,seal=seal)
    save(out/'assembly_audit.json',result)
    table(out/'fixture_contact_audit.csv',contacts);table(out/'closure_insertion.csv',[dict(kind='CLOSURE',**r) for r in closing]+[dict(kind='LOCK_PIN',**r) for r in locking]);table(out/'core_to_dish_clearance.csv',[dict(core_surface_min_z_mm=floor,dish_plane_z_mm=0,clearance_mm=floor)]);table(out/'core_to_wall_clearance.csv',contacts)
    table(out/'assembly_clearances.csv',[dict(port_id=p['port_id'],radial_mm=p['radial_clearance_mm'],slide_mm=p['slide_clearance_mm'],axial_float_plus_minus_mm=p['axial_float_mm']) for p in ports]);table(out/'port_fixture_design.csv',ports)
    table(out/'port_constraint_scheme.csv',[dict(port_id=p['port_id'],role=p['role'],radial_directions=['Z','wall_tangent'],axial_float_mm=p['axial_float_mm'],axial_stop=p['role']=='MASTER_DATUM',closure_lock='REMOVABLE_EXTERNAL_PIN') for p in ports]);table(out/'seal_geometry.csv',[dict(feature='BOTTOM_PERIMETER',closed=True,external=True),*[dict(feature='PORT_'+p['port_id'],closed=True,external=True) for p in ports]])
    return result
