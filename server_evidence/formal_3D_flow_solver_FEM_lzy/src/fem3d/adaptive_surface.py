"""Gmsh cap triangulation in rim-scaled coordinates; persistent 3D rim untouched."""
import time
import numpy as np
from .adaptive_port import rim_distance,size_function,surface_quality
from .cap_remesh import project,lift,cross2,check_wall
from .planar_port import validate_port
from .mesh_input import surface_topology


def triangulate_polygon(rim_xy,h_rim,H,policy):
    import gmsh
    scale=h_rim;xy=np.asarray(rim_xy)/scale
    gmsh.initialize()
    try:
        gmsh.model.add('port')
        options={'General.NumThreads':1,'General.Verbosity':1,'Mesh.MaxNumThreads1D':1,'Mesh.MaxNumThreads2D':1,
          'Mesh.Algorithm':6,'Mesh.ElementOrder':1,'Mesh.MeshSizeFromPoints':0,'Mesh.MeshSizeFromCurvature':0,
          'Mesh.MeshSizeExtendFromBoundary':0,'Mesh.Smoothing':0,'Mesh.MeshSizeMin':min(1.,H/scale),'Mesh.MeshSizeMax':max(1.,H/scale),
          'Geometry.Tolerance':policy['surface_meshing']['normalized_geometry_tolerance'],'Mesh.RandomSeed':policy['surface_meshing']['random_seed']}
        for k,v in options.items():gmsh.option.setNumber(k,v)
        cad=[gmsh.model.geo.addPoint(float(x),float(y),0,1.) for x,y in xy]
        lines=[gmsh.model.geo.addLine(cad[i],cad[(i+1)%len(cad)]) for i in range(len(cad))]
        face=gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(lines)]);gmsh.model.geo.synchronize()
        for line in lines:gmsh.model.mesh.setTransfiniteCurve(line,2)
        def sizing(dim,tag,x,y,z,lc):
            distance=rim_distance([[x,y]],xy)[0]
            return float(size_function(distance,1.,H/scale,policy['search']['grading_slope']))
        gmsh.model.mesh.setSizeCallback(sizing);gmsh.model.mesh.generate(2)
        node_to_local={};boundary_nodes=[]
        for i,tag in enumerate(cad):
            nodes,coords,_=gmsh.model.mesh.getNodes(0,tag)
            if len(nodes)!=1:raise ValueError('Rim CAD point node mismatch')
            node_to_local[int(nodes[0])]=i;boundary_nodes.append(int(nodes[0]))
        for i,line in enumerate(lines):
            if len(gmsh.model.mesh.getNodes(1,line,includeBoundary=False)[0]):raise ValueError('Rim split')
            types,_,nodes=gmsh.model.mesh.getElements(1,line)
            if list(types)!=[1] or set(map(int,nodes[0]))!={boundary_nodes[i],boundary_nodes[(i+1)%len(cad)]}:raise ValueError('Rim connectivity changed')
        nodes,coords,_=gmsh.model.mesh.getNodes();coords=coords.reshape(-1,3);result=list(np.asarray(rim_xy))
        for node,coord in sorted(zip(nodes,coords),key=lambda v:v[0]):
            if int(node) in node_to_local:continue
            node_to_local[int(node)]=len(result);result.append(coord[:2]*scale)
        types,_,nodes=gmsh.model.mesh.getElements(2,face)
        if list(types)!=[2]:raise ValueError('Non-triangular surface element')
        tri=np.array([node_to_local[int(i)] for i in nodes[0]],dtype=np.int64).reshape(-1,3)
        result=np.asarray(result);t=result[tri];flipped=cross2(t[:,1]-t[:,0],t[:,2]-t[:,0])<0;tri[flipped]=tri[flipped][:,[0,2,1]]
        return result,tri,{'gmsh_version':gmsh.__version__,'options':options,'coordinate_normalization_m':scale,'rim_geometry_reconstructed':False}
    finally:
        gmsh.finalize()


def vascular_trial(source,port,H,policy):
    start=time.perf_counter();ids=np.asarray(port['rim_vertex_ids']);rim=source['points_m'][ids]
    edge=source['points_m'][np.asarray(port['rim_edge_ids'])];h_rim=float(np.median(np.linalg.norm(edge[:,1]-edge[:,0],axis=1)))
    xy,_=project(rim,port['plane_origin_m'],np.array(port['basis']))
    local_xy,local_tri,metadata=triangulate_polygon(xy,h_rim,H,policy)
    interior=lift(local_xy[len(ids):],port['plane_origin_m'],np.array(port['basis']))
    n=len(source['points_m']);points=np.vstack([source['points_m'],interior]);index=np.r_[ids,np.arange(n,len(points))]
    cap=index[local_tri];keep=source['facet_tags']!=port['entity_id']
    combined={'points_m':points,'triangles':np.vstack([source['triangles'][keep],cap]),'facet_tags':np.r_[source['facet_tags'][keep],np.full(len(cap),port['entity_id'],dtype=np.int32)]}
    # Geometry precedes all quality decisions, including a full closed-surface audit.
    try:
        geometry=validate_port(source['points_m'],points,cap,port,np.arange(n,len(points)),policy,origin=port['plane_origin_m'],normal=port['outward_normal'],basis=port['basis'])
        wall=check_wall(source,combined);topology=surface_topology(points,combined['triangles'])
    except ValueError as exc:
        return {'points_m':points,'cap_triangles':cap}, {'geometry_status':'FAIL','geometry_error':str(exc),'quality':None,'H':H,'h_rim':h_rim,'R_eq':float(np.sqrt(port['formal_projected_area_m2']/np.pi)),'triangle_count':len(cap),'runtime_s':time.perf_counter()-start,'metadata':metadata}
    quality=surface_quality(points,cap,policy)
    data={'points_m':points,'cap_triangles':cap}
    return data,{'geometry_status':'PASS','geometry':geometry,'wall':wall,'combined_surface_topology':topology,'quality':quality,'H':H,'h_rim':h_rim,
       'R_eq':float(np.sqrt(port['formal_projected_area_m2']/np.pi)),'triangle_count':len(cap),'Steiner_count':len(interior),'runtime_s':time.perf_counter()-start,'metadata':metadata}


def combine_ports(source,contract,port_meshes,policy):
    points=source['points_m'].copy();n=len(points);tris=[source['triangles'][source['facet_tags']==1]];tags=[source['facet_tags'][source['facet_tags']==1]];results={}
    for name,port in contract['ports'].items():
        d=port_meshes[name];new=d['points_m'][n:];mapping=np.r_[np.arange(n),np.arange(len(points),len(points)+len(new))]
        tri=mapping[d['cap_triangles']];points=np.vstack([points,new]);tris.append(tri);tags.append(np.full(len(tri),port['entity_id'],np.int32))
        interior=np.setdiff1d(np.unique(tri),port['rim_vertex_ids'])
        geom=validate_port(source['points_m'],points,tri,port,interior,policy,origin=port['plane_origin_m'],normal=port['outward_normal'],basis=port['basis'])
        quality=surface_quality(points,tri,policy)
        results[name]={'geometry':geom,'quality':quality,'triangle_count':len(tri),'Steiner_count':len(new)}
    surface={'points_m':points,'triangles':np.vstack(tris),'facet_tags':np.concatenate(tags)}
    wall=check_wall(source,surface);topology=surface_topology(points,surface['triangles'])
    return surface,{'status':'PASS' if all(p['quality']['status']=='PASS' for p in results.values()) else 'FAIL','wall':wall,'topology':topology,'ports':results,'total_cap_triangles':sum(p['triangle_count'] for p in results.values()),'cap_count_is_hard_gate':False}
