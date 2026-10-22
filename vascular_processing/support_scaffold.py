"""Independent experimental pillar/short-bridge alternative (paper idea only)."""
import numpy as np
import trimesh
from .support_tree_planner import strut,illegal_collision

def bridge(a,b,diameter,max_span):
    length=float(np.linalg.norm(np.asarray(a)-b))
    if length>max_span:raise ValueError('BRIDGE_SPAN_EXCEEDS_ENGINEERING_LIMIT')
    return strut(a,b,diameter),dict(length_mm=length,limit_mm=max_span,physical_print_validated=False,status='EXPERIMENTAL_BRIDGE')

def build(contacts,vascular,box,inner,profile,diameter=2.5,max_span=6.):
    pieces=[];edges=[];rows=[];pairs=[]
    for i,c in enumerate(contacts):
        top=np.array(c['xyz'])-[0,0,profile['support_top_gap']];shoulder=top-[0,0,1.5];base=shoulder.copy()
        inside=np.all(base[:2]>inner[0,:2]) and np.all(base[:2]<inner[1,:2]);base[2]=inner[0,2]+.2 if inside else 0
        m=trimesh.boolean.union([strut(base,shoulder,diameter,1.2),strut(shoulder,top,1.2,c.get('tip_diameter',.6))],engine='manifold')
        q=illegal_collision(m,vascular,box);rows.append(dict(contact_id=c['contact_id'],**q))
        pieces.append(m);edges.append((base,shoulder,top))
    for i in range(0,len(edges)-1,2):
        a=edges[i][1].copy();b=edges[i+1][1].copy();a[2]=b[2]=min(a[2],b[2])-1
        # Different anatomical branches can share an attachment position.
        # Their pillars already coincide; a zero-length connector is neither
        # a bridge nor an error in the entire research parameter sweep.
        if 1e-8 < np.linalg.norm(a-b) <= max_span:
            m,r=bridge(a,b,1.,max_span);q=illegal_collision(m,vascular,box);r.update(q);pairs.append(r)
            if q['passed']:pieces.append(m)
    mesh=trimesh.boolean.union(pieces,engine='manifold') if pieces else trimesh.Trimesh()
    return mesh,dict(pillars=rows,bridges=pairs,diameter_mm=diameter,max_bridge_span_mm=max_span,production=False,status='RESEARCH_CANDIDATE',physical_print_validated=False)

def self_printability(mesh,profile,bridge_limit=0.):
    from .support_detection import layers
    result=layers(mesh,profile,.01);uncovered=sum(r['unsupported_area_mm2'] for r in result['rows'])
    # Disconnected first layer on the mold floor is intentionally unresolved
    # until slicing together with the frozen model. Never label it printable
    # merely because the support mesh is watertight.
    return dict(unsupported_layer_region_count=len(result['rows']),unsupported_layer_area_sum_mm2=uncovered,passed=not result['rows'],
        bridge_limit_mm=bridge_limit,interpretation='Layer self-support only; any bridge or elevated base still requires actual combined slicing')
