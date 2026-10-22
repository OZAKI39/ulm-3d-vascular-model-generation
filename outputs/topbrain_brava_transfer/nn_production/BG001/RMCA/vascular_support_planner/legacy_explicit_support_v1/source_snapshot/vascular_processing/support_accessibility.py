"""Top access of actual support components using the existing FCL capsule probe."""
import numpy as np
import trimesh
from .support_access_qc import capsule,line_visible

def component_access(components,vascular,box,outer,probe_diameter=3.):
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular',vascular);manager.add_object('box',box)
    rows=[];top=float(np.asarray(outer)[1,2]);r=probe_diameter/2
    for i,component in enumerate(components):
        # Grasp points on four sides of the component, including low trunks.
        bounds=component.bounds;targets=[]
        actual=np.asarray(component.metadata.get('support_occupied_points',component.vertices))
        for z in np.quantile(actual[:,2],[.2,.5,.8]):
            band=actual[np.abs(actual[:,2]-z)<=max(.3,(bounds[1,2]-bounds[0,2])*.1)]
            if not len(band):band=actual[[int(np.argmin(abs(actual[:,2]-z)))]]
            for axis in [0,1]:
                for side in [0,1]:
                    index=np.argmin(band[:,axis]) if side==0 else np.argmax(band[:,axis])
                    p=band[index].copy();p[axis]+=(-1 if side==0 else 1)*(r+.25);targets.append(p)
        visible=False;reachable=False;accepted=None
        for p in targets:
            for offset in [[0,0],[3,0],[-3,0],[0,3],[0,-3]]:
                entry=np.r_[p[:2]+offset,top+r+2];ray=line_visible(box,vascular,entry,p);visible|=ray
                if ray and not manager.in_collision_single(capsule(entry,p,r,12)):
                    reachable=True;accepted=dict(start=entry.tolist(),end=p.tolist());break
            if reachable:break
        rows.append(dict(component_id=i,top_reachable=visible,probe_3mm_reachable=reachable,trapped_support_risk=not reachable,
            enclosed_and_invisible=not visible,status='TOP_ACCESS_FOUND' if reachable else 'TRAPPED_SUPPORT_RISK',probe_path=accepted,
            limitation='Sampled straight/small-angle 3 mm approach; grasp access does not certify intact removal'))
    return dict(rows=rows,trapped_support_count=sum(r['trapped_support_risk'] for r in rows),top_accessible_fraction=float(np.mean([r['probe_3mm_reachable'] for r in rows])) if rows else None,
        method='FCL capsule against frozen vascular + box; grasp targets adjacent to actual support vertices/occupied extrusion cells, never empty bounding-box corners')
