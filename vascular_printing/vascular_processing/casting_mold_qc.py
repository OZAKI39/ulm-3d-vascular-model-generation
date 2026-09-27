"""Solid preservation, wall passages and topology-aware triangle-surface QC."""
import itertools
import numpy as np
import networkx as nx
from scipy.spatial import cKDTree
import trimesh
from . import sacrificial_fixture as old
from .sacrificial_print_frame import box_mesh,port_stem_bounds
from .print_frame_qc import intersection,volume,interface_geometry_checks


def intact_box_qc(box,cfg):
    mesh=box['mesh'];inner=box['inner'];outer=box['outer'];rows=[]
    for name,bounds in box['members'].items():
        expected=box_mesh(bounds);missing=trimesh.boolean.difference([expected,mesh],engine='manifold')
        rows.append(dict(wall=name,missing_material_mm3=volume(missing),passed=volume(missing)<=cfg['geometry']['volume_tolerance_mm3']))
    cavity=box_mesh(inner);cavity_overlap=volume(intersection(mesh,cavity))
    # Rays start within the cavity, traverse its top and stop outside. They
    # distinguish an open material boundary from a falsely sealed roof.
    xy=np.array(list(itertools.product(np.linspace(inner[0,0]+.1,inner[1,0]-.1,9),np.linspace(inner[0,1]+.1,inner[1,1]-.1,9))))
    origins=np.column_stack([xy,np.full(len(xy),inner[1,2]-.5)])
    hits,ids,_=mesh.ray.intersects_location(origins,np.tile([0.,0,1.],(len(origins),1)),multiple_hits=True)
    blocks=sum(h[2]>=inner[1,2]-.5 for h in hits)
    dimension=inner[1]-inner[0]
    return dict(walls=rows,four_side_walls_closed=all(r['passed'] for r in rows if r['wall']!='BOTTOM'),
        bottom_closed=rows[-1]['passed'],top_fully_open=blocks==0 and cavity_overlap<=cfg['geometry']['volume_tolerance_mm3'],
        upward_test_ray_count=len(origins),blocked_top_ray_count=int(blocks),cavity_box_overlap_mm3=cavity_overlap,
        top_opening_width_mm=float(dimension[0]),top_opening_length_mm=float(dimension[1]),top_opening_area_mm2=float(dimension[0]*dimension[1]),
        inner_dimensions_mm=dimension,outer_dimensions_mm=outer[1]-outer[0],wall_thickness_mm=box['wall_thickness_mm'],
        bottom_thickness_mm=box['bottom_thickness_mm'],dimensions_policy=box['dimensions_policy'],
        no_assembly_clearance_holes=all(r['passed'] for r in rows),mesh=old.mesh_qc(mesh,cfg))


def port_overlaps(inputs,box,cfg):
    contacts=intersection(inputs['core'],box['mesh']);rows=[];parts={}
    for p in inputs['ports']:
        contact=intersection(contacts,box_mesh(port_stem_bounds(p,cfg['geometry']['coordinate_tolerance_mm'])))
        v=volume(contact);direction=p['direction'];projections=contact.vertices@direction if len(contact.vertices) else []
        span=float(np.ptp(projections)) if len(projections) else 0.
        axis,side=old.face_axis(p['face']);entry=p['target'].copy();exit=entry.copy()
        entry[axis]=box['inner'][side,axis];exit[axis]=box['outer'][side,axis]
        okay=v>0 and span>=cfg['box']['minimum_overlap_fraction']*box['wall_thickness_mm']
        rows.append(dict(port_id=p['port_id'],wall=p['face'],port_radius_mm=p['radius_mm'],
            wall_thickness_mm=box['wall_thickness_mm'],intersection_volume_mm3=v,intersection_axial_span_mm=span,
            wall_entry_point=entry.tolist(),wall_exit_point=exit.tolist(),status='PASS' if okay else 'PORT_WALL_OVERLAP_INSUFFICIENT_'+p['port_id']))
        parts[p['port_id']]=contact
    return rows,parts


def union_audit(core,box,result,cfg):
    contact=intersection(core,box)
    core_loss=volume(trimesh.boolean.difference([core,result],engine='manifold'))
    box_loss=volume(trimesh.boolean.difference([box,result],engine='manifold'))
    result_qc=old.mesh_qc(result,cfg);failures=[]
    if core_loss>cfg['geometry']['volume_tolerance_mm3']:failures.append('VASCULAR_CORE_DAMAGED_BY_BOX_UNION')
    if box_loss>cfg['geometry']['volume_tolerance_mm3']:failures.append('BOX_MATERIAL_DAMAGED_BY_UNION')
    if not result_qc['passed']:failures.append('FINAL_MOLD_MESH_QC_FAILED')
    return dict(core_volume_before_mm3=volume(core),box_volume_before_mm3=volume(box),intersection_volume_mm3=volume(contact),
        union_volume_mm3=volume(result),volume_balance_error_mm3=abs(volume(result)-(volume(core)+volume(box)-volume(contact))),
        core_material_loss_mm3=core_loss,box_material_loss_mm3=box_loss,**result_qc,euler=int(result.euler_number),failures=failures)


def closest_surface(mesh,points):
    import vtk
    from .sacrificial_fixture_review import polydata
    locator=vtk.vtkStaticCellLocator();locator.SetDataSet(polydata(mesh));locator.BuildLocator()
    near=[];distance=[];triangles=[];p=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);squared=vtk.reference(0.)
    for query in points:
        locator.FindClosestPoint(query,p,cell,sub,squared)
        near.append(p.copy());distance.append(np.sqrt(float(squared)));triangles.append(int(cell))
    return np.array(near),np.array(distance),np.array(triangles,dtype=int)


def companion_samples(inputs,cfg):
    data=inputs['companion'];graph=data['graph'];transform=data['transform'];samples=[];labels=[];origins=[]
    for a,b in graph.edges:
        pa,pb=graph.nodes[a]['coords'][:3],graph.nodes[b]['coords'][:3];length=np.linalg.norm(pb-pa)
        for t in np.linspace(0,1,max(2,int(np.ceil(length/cfg['geometry']['centerline_label_sample_step_mm']))+1)):
            samples.append(transform[:3,:3]@(pa*(1-t)+pb*t)+transform[:3,3]);labels.append(data['branches'][a,b]);origins.append((a,b,t,length))
    return np.array(samples),np.array(labels),origins


def partition_vascular(final,inputs,box,cfg):
    """Attribute actual union surface using immutable geometry and SWC identity.

    Proximity only assigns identities. Every reported separation uses FCL
    triangle surfaces, never sampled centerline/vertex distances.
    """
    centers=final.triangles_center;_,d,_=closest_surface(inputs['core'],centers)
    vascular_ids=np.flatnonzero(d<cfg['geometry']['coordinate_tolerance_mm'])
    points=centers[vascular_ids]
    meshes=[inputs['central']]+[p['reference_mesh'] for p in inputs['ports']]
    distances=np.column_stack([closest_surface(m,points)[1] for m in meshes])
    owner=np.argmin(distances+np.arange(len(meshes))[None,:]*1e-9,axis=1)
    samples,branch_labels,origins=companion_samples(inputs,cfg);tree=cKDTree(samples)
    labels=branch_labels[tree.query(points)[1]].astype('U40')
    for j,p in enumerate(inputs['ports'],1):labels[owner==j]=p['port_id']
    interior=np.all((final.triangles[vascular_ids]>=box['inner'][0]-1e-5)&(final.triangles[vascular_ids]<=box['inner'][1]+1e-5),axis=(1,2))
    groups={name:vascular_ids[(labels==name)&interior] for name in sorted(set(labels))}
    inner_walls={}
    for axis,side,name in [(0,0,'-X'),(0,1,'+X'),(1,0,'-Y'),(1,1,'+Y'),(2,0,'BOTTOM')]:
        normal=np.zeros(3);normal[axis]=1 if side==0 else -1
        ids=np.flatnonzero((abs(centers[:,axis]-box['inner'][side,axis])<1e-5)&(final.face_normals@normal>.99))
        inner_walls['BOX_INNER_'+name]=ids
    return dict(face_ids=vascular_ids,labels=labels,interior=interior,groups=groups,samples=samples,
        branch_labels=branch_labels,sample_origins=origins,sample_tree=tree,
        inner_wall_groups=inner_walls,
        method='Frozen component-surface attribution + companion SWC labels; actual delivered union triangles')


def fcl_distance(a,b):
    manager=trimesh.collision.CollisionManager();manager.add_object('first',a)
    value,details=manager.min_distance_single(b,return_data=True)
    return float(max(0,value)),np.array(details.point('first')),np.array(details.point('__external'))


def classify_gap(distance):
    return 'COMFORTABLE' if distance>=3 else 'MODERATE' if distance>=2 else 'THIN_LOCAL' if distance>=1.5 else 'VERY_THIN'


def topology_exemptions(inputs,partition):
    graph=inputs['companion']['graph'];branches=inputs['companion']['branches'];branch_nodes={}
    for (a,b),label in branches.items():branch_nodes.setdefault(label,set()).update((a,b))
    weighted=graph.to_undirected()
    for a,b in weighted.edges:weighted[a][b]['weight']=float(np.linalg.norm(graph.nodes[a]['coords'][:3]-graph.nodes[b]['coords'][:3]))
    return branch_nodes,weighted


def trim_near_node(patch,node,partition,weighted,cfg):
    """Clip analysis surfaces at the geodesic exemption boundary.

    Removing whole triangles would leave a false zero at a junction when a
    coarse synthetic triangle spans both sides of the exemption boundary.
    VTK clips only this diagnostic patch; the production STL is unchanged.
    """
    if not len(patch.faces):return patch
    from .sacrificial_fixture_review import polydata
    graph_distance=nx.single_source_dijkstra_path_length(weighted,node,weight='weight')
    distances=np.array([min(graph_distance[a]+t*length,graph_distance[b]+(1-t)*length) for a,b,t,length in partition['sample_origins']])
    values=distances[partition['sample_tree'].query(patch.vertices)[1]]
    data=polydata(patch);data.point_data['junction_arclength_mm']=values
    clipped=data.clip_scalar(scalars='junction_arclength_mm',value=cfg['pdms']['junction_exemption_mm'],invert=False).extract_surface(algorithm='dataset_surface').triangulate()
    if not clipped.n_cells:return trimesh.Trimesh(vertices=np.empty((0,3)),faces=np.empty((0,3),int),process=False)
    return trimesh.Trimesh(clipped.points,clipped.faces.reshape(-1,4)[:,1:],process=False)


def ligament_qc(final,inputs,partition,cfg):
    groups=partition['groups'];branches,graph=topology_exemptions(inputs,partition);ports={p['port_id']:p for p in inputs['ports']};rows=[]
    def sub(ids):return final.submesh([ids],append=True,repair=False)
    for first,second in itertools.combinations(sorted(groups),2):
        a=sub(groups[first]);b=sub(groups[second]);excluded='NONE'
        if first in branches and second in branches:
            shared=branches[first]&branches[second]
            if shared:
                excluded='LOCAL_TRUE_BIFURCATION'
                for node in shared:
                    a=trim_near_node(a,node,partition,graph,cfg);b=trim_near_node(b,node,partition,graph,cfg)
        elif (first in ports) != (second in ports):
            pname=first if first in ports else second;nname=second if first in ports else first
            if ports[pname]['parent_branch']==nname:
                excluded='LOCAL_PORT_PARENT_ATTACHMENT'
                if first==nname:a=trim_near_node(a,ports[pname]['swc_id'],partition,graph,cfg)
                else:b=trim_near_node(b,ports[pname]['swc_id'],partition,graph,cfg)
        if not len(a.faces) or not len(b.faces):continue
        distance,p,q=fcl_distance(a,b)
        rows.append(dict(first=first,second=second,distance_mm=distance,geometry_grade=classify_gap(distance),
            first_point_mm=p.tolist(),second_point_mm=q.tolist(),exclusion=excluded,method='FCL_TRIANGLE_SURFACE'))
    rows.sort(key=lambda r:(r['distance_mm'],r['first'],r['second']))
    o3=min((r for r in rows if 'O3' in (r['first'],r['second'])),key=lambda r:r['distance_mm'],default=None)
    baseline=inputs.get('summary',{}).get('per_port',{}).get('O3',{}).get('clearance_mm')
    passed=o3 is not None and baseline is not None and o3['distance_mm']>=cfg['pdms']['o3_minimum_mm'] and abs(o3['distance_mm']-baseline)<=cfg['pdms']['o3_maximum_deviation_mm']
    return dict(minimum_ligament_mm=rows[0]['distance_mm'],minimum_pair=rows[0],top10=rows[:10],all_pairs=rows,
        o3_gap_mm=o3['distance_mm'] if o3 else None,o3_pair=o3,o3_accepted_gap_mm=baseline,o3_deviation_mm=abs(o3['distance_mm']-baseline) if o3 and baseline is not None else None,
        o3_status=('PASS' if passed else 'O3_PDMS_GAP_REGRESSION') if 'O3' in ports else 'NOT_APPLICABLE_SYNTHETIC',same_branch_pairs_excluded=True,
        junction_exemption_mm=cfg['pdms']['junction_exemption_mm'],exclusion_method='Only local true junction/port-parent regions are excluded; distal portions remain checked.')


def wall_clearance(final,partition,inputs,box,cfg):
    rows=[];ports={p['port_id']:p for p in inputs['ports']}
    for name,ids in partition['groups'].items():
        patch=final.submesh([ids],append=True,repair=False)
        if not len(patch.faces):continue
        for wall,bounds in box['members'].items():
            current=patch;exempt=False
            if name in ports and ports[name]['face']==wall:
                p=ports[name];origin=p['target']-p['direction']*cfg['pdms']['wall_crossing_exemption_mm']
                current=patch.slice_plane(origin,-p['direction'],cap=False);exempt=True
            if not len(current.faces):continue
            # For an interior vascular patch, the nearest surface of this wall
            # solid is its true inner face (bottom included).
            distance,a,b=fcl_distance(current,box_mesh(bounds))
            rows.append(dict(vascular_branch=name,wall=wall,distance_mm=distance,vascular_point_mm=a.tolist(),wall_point_mm=b.tolist(),
                legal_port_crossing_excluded=exempt,excluded_axial_band_mm=cfg['pdms']['wall_crossing_exemption_mm'] if exempt else 0.,
                method='FCL_TRIANGLE_SURFACE',interpretation='Distance outside the documented legal penetration band'))
    return sorted(rows,key=lambda r:(r['distance_mm'],r['vascular_branch'],r['wall']))
