"""CPU port of FreeFloating's support-point graph, not its renderer/topology.

Adapted from regrow1123/FreeFloating, supportpoint.cpp/h, commit
5462f823abf01661e9abeb85d9bc871e01d4ec39. Copyright (c) 2019 Seongje Jang.
MIT license: third_party/support_research/FreeFloating/LICENSE.
OpenGL fragment lists are replaced by trimesh/Embree ray intersections.
The original entry/exit graph and quantized layer rule are retained. Dijkstra
termination is made explicit on settled distances (rather than throwing from
Boost's discover visitor), and ground-level points do not require contacts.
"""
import heapq
import numpy as np

def machine_layer(z,height):return np.ceil((np.asarray(z)-height*.5)/height).astype(int)

def graph_from_intervals(cells,height,angle=np.pi/4):
    xyz=[];indices={}
    for cell,values in sorted(cells.items()):
        ids=[]
        for p in values:ids.append(len(xyz));xyz.append(p)
        indices[cell]=ids
    xyz=np.asarray(xyz).reshape(-1,3);edges=[[] for _ in xyz]
    for (i,j),ids in sorted(indices.items()):
        for k in range(0,len(ids)-1,2):
            target=ids[k];p=xyz[target];layer=machine_layer(p[2],height)
            edges[target].append((ids[k+1],0.))
            for di,dj in [(a,b) for a in [-1,0,1] for b in [-1,0,1] if a or b]:
                other=indices.get((i+di,j+dj),[])
                for q in range(0,len(other)-1,2):
                    entry,exit=xyz[other[q]],xyz[other[q+1]]
                    if machine_layer(entry[2],height)>layer:break
                    if machine_layer(exit[2],height)>=layer:
                        v=p-entry;theta=np.arccos(np.clip(v[2]/max(np.linalg.norm(v),1e-12),-1,1))
                        weight=float(np.clip((theta-angle)/(np.pi/2-angle),0,1))
                        edges[other[q]].append((target,weight));break
    return xyz,edges

def select_graph_points(xyz,edges,coverage,ground_z=0.):
    marked=np.zeros(len(xyz),bool);selected=[]
    for source in sorted(range(len(xyz)),key=lambda i:(xyz[i,2],i)):
        if marked[source]:continue
        if xyz[source,2]>ground_z+1e-5:selected.append(source)
        queue=[(0.,source)];dist={source:0.}
        while queue:
            d,i=heapq.heappop(queue)
            if d>coverage:break
            if d>dist[i]:continue
            marked[i]=True
            for j,w in edges[i]:
                nd=d+w
                if nd<=coverage and nd<dist.get(j,float('inf')):dist[j]=nd;heapq.heappush(queue,(nd,j))
    return xyz[selected],selected

def detect(mesh,layer_height,pixel=.35,effective_radius=1.):
    lo,hi=mesh.bounds;x=np.arange(lo[0]-pixel*.31,hi[0]+pixel,pixel);y=np.arange(lo[1]-pixel*.27,hi[1]+pixel,pixel)
    grid=np.stack(np.meshgrid(x,y,indexing='ij'),axis=-1).reshape(-1,2)
    origins=np.column_stack([grid,np.full(len(grid),lo[2]-1)])
    hits,rays,_=mesh.ray.intersects_location(origins,np.tile([0,0,1],(len(grid),1)),multiple_hits=True)
    groups={}
    for p,r in zip(hits,rays):groups.setdefault(int(r),[]).append(p)
    cells={};odd=0
    for r,points in sorted(groups.items()):
        points=np.array(sorted(points,key=lambda p:p[2]));points=points[np.r_[True,np.diff(points[:,2])>1e-5]]
        if len(points)%2:odd+=1;continue
        cells[(r//len(y),r%len(y))]=points
    xyz,edges=graph_from_intervals(cells,layer_height)
    points,ids=select_graph_points(xyz,edges,effective_radius/pixel)
    return points,dict(implementation='MIT_CPU_PORT_RAY_INTERVAL_GRAPH',nodes=len(xyz),edges=sum(map(len,edges)),odd_ray_columns_excluded=odd,
        ray_columns=len(grid),pixel_mm=pixel,effective_radius_mm=effective_radius,selected_count=len(ids),reference_binary_status='BUILD_UNAVAILABLE_MISSING_BOOST_OPENGL_DEPENDENCIES',
        differences=['CPU ray intervals replace OpenGL rasterizer','Deterministic ties','Coverage cutoff on settled Dijkstra distance'])

def local_minima(mesh,tolerance=1e-6):
    z=mesh.vertices[:,2];ids=[]
    for i,neighbors in enumerate(mesh.vertex_neighbors):
        other=z[neighbors]
        if len(other) and np.all(z[i]<=other+tolerance) and np.any(z[i]<other-tolerance):ids.append(i)
    return mesh.vertices[ids]
