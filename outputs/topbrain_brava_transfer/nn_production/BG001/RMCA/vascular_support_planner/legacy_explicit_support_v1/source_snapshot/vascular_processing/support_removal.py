"""Whole-component deterministic removal after neck cutting, with swept bounds."""
from __future__ import annotations
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation,Slerp

def cut_necks(mesh,roots):
    if not roots:return []
    z=max(r['estimated_break_location'][2] for r in roots)
    lo,hi=mesh.bounds;size=hi-lo+4;size[2]=z-lo[2]+2
    cutter=trimesh.creation.box(size);cutter.apply_translation([(lo[0]+hi[0])/2,(lo[1]+hi[1])/2,z-size[2]/2])
    broken=trimesh.boolean.difference([mesh,cutter],engine='manifold')
    return list(broken.split(only_watertight=False))

def swept_path(component,obstacles,poses,tolerance=.004,static_manager=None):
    """FCL distances plus a rigid-motion Lipschitz bound between poses.

    Endpoint-only tests miss thin obstacles. Every interval must have distance
    greater than its maximum vertex travel, or is subdivided; unresolved narrow
    clearance fails conservatively. This is not a robotics motion planner.
    """
    manager=static_manager
    if manager is None:
        manager=trimesh.collision.CollisionManager()
        for i,m in enumerate(obstacles):manager.add_object(str(i),m)
    moving=trimesh.collision.CollisionManager();moving.add_object('component',component)
    radius=float(np.linalg.norm(component.vertices,axis=1).max());cache={};minimum=np.inf;count=0
    def sample(m):
        nonlocal minimum,count
        key=m.tobytes()
        if key not in cache:
            moving.set_transform('component',m);hit=manager.in_collision_other(moving)
            dist=0. if hit else float(manager.min_distance_other(moving));cache[key]=(hit,dist);count+=1;minimum=min(minimum,dist)
        return cache[key]
    def interval(a,b,depth=0):
        hit_a,da=sample(a);hit_b,db=sample(b)
        if hit_a or hit_b:return False,'SWEPT_COLLISION'
        angle=Rotation.from_matrix(a[:3,:3].T@b[:3,:3]).magnitude();travel=float(np.linalg.norm(b[:3,3]-a[:3,3])+angle*radius)
        if min(da,db)>travel+.0001:return True,'CLEAR'
        if travel<tolerance or depth>=18:return False,'CLEARANCE_NOT_CERTIFIED'
        middle=np.eye(4);middle[:3,:3]=Slerp([0,1],Rotation.from_matrix(np.array([a[:3,:3],b[:3,:3]])))(.5).as_matrix();middle[:3,3]=(a[:3,3]+b[:3,3])/2
        passed,status=interval(a,middle,depth+1)
        return interval(middle,b,depth+1) if passed else (False,status)
    for a,b in zip(poses[:-1],poses[1:]):
        ok,status=interval(a,b)
        if not ok:return dict(passed=False,status=status,min_clearance_mm=float(minimum),distance_queries=count)
    return dict(passed=True,status='SWEPT_PATH_CLEAR',min_clearance_mm=float(minimum),distance_queries=count)

def removal_paths(component,top,angles=(10,20)):
    pivot=np.r_[component.centroid[:2],component.bounds[0,2]];local=component.copy();local.apply_translation(-pivot)
    initial=np.eye(4);initial[:3,3]=pivot;rise=float(top-component.bounds[0,2]+5)
    up=initial.copy();up[2,3]+=rise
    paths=[('VERTICAL_PLUS_Z',[initial,up])]
    for axis in ['x','y']:
        for degrees in angles:
            for sign in [-1,1]:
                tilted=initial.copy();tilted[:3,:3]=Rotation.from_euler(axis,sign*degrees,degrees=True).as_matrix()
                final=tilted.copy();final[2,3]+=rise
                paths.append((f'TILT_{axis}_{sign*degrees:+g}_THEN_UP',[initial,tilted,final]))
                # Tilting only the top does not free the lower trunk's swept
                # projection: it can still hit the vessel later on ascent.
                # Release the whole component sideways in the lean direction
                # before lifting. The same FCL swept test checks every leg.
                direction=tilted[:3,:3]@np.array([0.,0.,1.]);xy=direction[:2]/np.linalg.norm(direction[:2])
                released=tilted.copy();released[:2,3]+=xy*2.;raised=released.copy();raised[2,3]+=rise
                paths.append((f'TILT_{axis}_{sign*degrees:+g}_LATERAL_2_THEN_UP',[initial,tilted,released,raised]))
    # Limited lateral release remains deterministic and explicit. First tilt
    # gently away from the tip, then translate at most 2 mm and move upward.
    for axis in [0,1]:
        for sign in [-1,1]:
            released=initial.copy();released[axis,3]+=sign*2.;up=released.copy();up[2,3]+=rise
            paths.append((f'LATERAL_{axis}_{sign*2:+g}_THEN_UP',[initial,released,up]))
    # A diagonal release can avoid a crossing branch when either pure X or Y
    # catches it. Keep the same bounded 2 mm distance, without random search.
    for axis in ['x','y']:
        for degrees in angles:
            for sign in [-1,1]:
                tilted=initial.copy();tilted[:3,:3]=Rotation.from_euler(axis,sign*degrees,degrees=True).as_matrix()
                for direction in [45,135,225,315]:
                    released=tilted.copy();released[:2,3]+=2*np.array([np.cos(np.radians(direction)),np.sin(np.radians(direction))])
                    raised=released.copy();raised[2,3]+=rise
                    paths.append((f'TILT_{axis}_{sign*degrees:+g}_DIAGONAL_{direction}_2_THEN_UP',[initial,tilted,released,raised]))
    return local,pivot,paths

def evaluate(mesh,roots,tips,vascular,box):
    components=cut_necks(mesh,roots);rows=[];manager=trimesh.collision.CollisionManager();manager.add_object('vascular',vascular);manager.add_object('box',box)
    for i,part in enumerate(components):
        local,pivot,paths=removal_paths(part,box.bounds[1,2]);attempts=[];accepted=None
        for name,poses in paths:
            result=swept_path(local,[vascular,box],poses,static_manager=manager);attempts.append(dict(path=name,**result))
            if result['passed']:accepted=dict(name=name,transforms=[m.tolist() for m in poses]);break
        lo,hi=part.bounds
        tip_ids=[p['tip_id'] for p in tips if p['root_connected'] and np.all(np.array(p['xyz'])>=lo-.01) and np.all(np.array(p['xyz'])<=hi+.01)]
        root_ids=[r['root_id'] for r in roots if np.all(np.array(r['xyz'])[:2]>=lo[:2]-.1) and np.all(np.array(r['xyz'])[:2]<=hi[:2]+.1)]
        rows.append(dict(component_id=i,tips=tip_ids,roots=root_ids,volume_mm3=float(part.volume),bbox_mm=part.bounds.tolist(),root_neck=[r for r in roots if r['root_id'] in root_ids],
            vertical_removable=attempts[0]['passed'],tilt_path=accepted,min_clearance_mm=attempts[-1]['min_clearance_mm'] if accepted else min(a['min_clearance_mm'] for a in attempts),
            status='REMOVAL_PATH_FOUND' if accepted else 'TRAPPED_SUPPORT_COMPONENT',attempts=attempts,pivot_mm=pivot.tolist(),
            physical_removal_validated=False,other_support_components_ignored=True))
    return components,dict(rows=rows,trapped_support_count=sum(r['status']=='TRAPPED_SUPPORT_COMPONENT' for r in rows),
        component_count=len(rows),method='FCL whole-component distance + bounded swept motion, after explicit neck removal',physical_removal_validated=False)
