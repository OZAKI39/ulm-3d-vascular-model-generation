"""Permanent validation geometries and evidence helpers, no production choices."""
from pathlib import Path
from dataclasses import asdict
import csv,json
import numpy as np
from .wall_geometry import WallGeometry
from .particle_shapes import Sphere,Ellipsoid,Capsule
from .rbc_distribution import sample_rbc_geometries,quantile_indices,stratified_indices
from .rbc import RBCGeometry
from .convex_triangle import triangle_gap
from .rbc_orientation import quaternion_from_short_axis,rotation_matrix

P2_COMMIT='c704e08d39f6134da300fc91cc93c8561314a17b'
P3_BRANCH='dev/particle-3-wall-contact-20260920'
P3_SELECTION_SEED=2026092128
P3_TUBE_RADII_M=np.linspace(.6,5.,9)*1e-6


def json_safe(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,Path):return str(value)
    raise TypeError(type(value))


def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False,default=json_safe)+'\n')


def write_rows(path,rows):
    rows=list(rows);path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=keys,lineterminator='\n');writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v,ensure_ascii=False,allow_nan=False,default=json_safe) if isinstance(v,(dict,list,tuple,np.ndarray)) else v for k,v in row.items()})


def population_and_geometries():
    p=sample_rbc_geometries(100000,2026092002)
    return p,[RBCGeometry.from_population(p,i) for i in quantile_indices(p.samples)]


def plane_triangle():return np.array([[-20,-20,0],[0,20,0],[20,-20,0]],dtype=float)*1e-6


def cylinder_wall(radius_m,sides=48,half_length_m=50e-6):
    theta=np.arange(sides)*2*np.pi/sides
    xy=np.column_stack((radius_m*np.cos(theta),radius_m*np.sin(theta)))
    lower=np.column_stack((xy,np.full(sides,-half_length_m)))
    upper=np.column_stack((xy,np.full(sides,half_length_m)))
    triangles=[]
    for i in range(sides):
        j=(i+1)%sides
        triangles.extend([[lower[i],lower[j],upper[j]],[lower[i],upper[j],upper[i]]])
    return WallGeometry(triangles,provenance=dict(role='VALIDATION_GEOMETRY_PARAMETERS',tube_circumradius_m=radius_m,
                        tube_apothem_m=radius_m*np.cos(np.pi/sides),sides=sides,half_length_m=half_length_m,open_endcaps=True))


def representative_wall_ids(wall,count=18):
    centers=wall.triangles.mean(axis=1);extent=np.ptp(centers,axis=0)
    normalized=(centers-centers.mean(axis=0))/np.maximum(extent,np.finfo(float).tiny)
    features=np.column_stack((normalized,wall.normal_in,.15*np.log(wall.areas_m2/np.median(wall.areas_m2))))
    chosen=[int(np.argmin(wall.areas_m2)),int(np.argmax(wall.areas_m2))]
    distance=np.min([np.sum((features-features[i])**2,axis=1) for i in chosen],axis=0)
    while len(chosen)<count:
        i=int(np.argmax(distance));chosen.append(i)
        distance=np.minimum(distance,np.sum((features-features[i])**2,axis=1))
    return chosen


def normal_audit(wall,field,ids):
    node_cells=[[] for _ in range(len(field.points_m))]
    for i,cell in enumerate(field.tetra):
        for node in cell:node_cells[int(node)].append(i)
    records=[]
    for i in ids:
        gids=wall.global_node_ids[i]
        cells=set(node_cells[int(gids[0])]).intersection(node_cells[int(gids[1])],node_cells[int(gids[2])])
        if len(cells)!=1:raise ValueError('WALL triangle must have exactly one owning lumen tetra')
        owner=cells.pop();center=field.points_m[field.tetra[owner]].mean(axis=0)
        tri=wall.triangles[i];centroid=tri.mean(axis=0)
        signed=float((center-centroid)@wall.normal_in[i])
        if signed<=0:raise ValueError('Frozen winding does not point inward to owning tetra')
        records.append(dict(triangle_id=i,vertices_m=tri,area_m2=wall.areas_m2[i],normal_out=wall.normal_out[i],
                            normal_in=wall.normal_in[i],owning_tetra=owner,owner_center_m=center,owner_inward_distance_m=signed))
    return records


def geometry_validation_records():
    records=[];tri=np.array([[0,0,0],[0,10,0],[10,0,0]],dtype=float)*1e-6
    for label,xy in [('FACE',[2,2]),('EDGE',[4,-.3]),('VERTEX',[-.3,-.4])]:
        for z in [.2,.6,1.,2.]:
            center=np.r_[xy,z]*1e-6;shape=Sphere(center,1e-6);g=triangle_gap(shape,tri)
            expected=np.linalg.norm(center-np.r_[np.maximum(xy,0),0]*1e-6)-1e-6
            records.append(dict(particle_type='MB',shape_mode=shape.mode,center_m=center,radius_m=1e-6,triangle_m=tri,
                feature=g.wall_feature,requested_feature=label,gap_m=g.gap_m,analytic_gap_m=expected,error_m=abs(g.gap_m-expected),
                wall_point_m=g.wall_point_m,particle_point_m=g.particle_point_m,normal_inward=g.normal_inward))
    ellipsoids=[];capsules=[];tri=plane_triangle()
    for ratio in [.15,.3,.6]:
        for angle in np.linspace(0,np.pi/2,13):
            r=np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]])
            for height in [.3e-6,1.5e-6,3e-6]:
                shape=Ellipsoid([0,0,height],[2e-6,2e-6,ratio*2e-6],r);g=triangle_gap(shape,tri)
                exact=height-2e-6*np.sqrt(np.sin(angle)**2+ratio**2*np.cos(angle)**2)
                ellipsoids.append(dict(shape_mode=shape.mode,ratio=ratio,angle_rad=angle,center_m=shape.center_m,
                    axes_m=shape.axes_m,rotation=r,gap_m=g.gap_m,analytic_gap_m=exact,error_m=abs(g.gap_m-exact),
                    triangle_m=tri,wall_point_m=g.wall_point_m,particle_point_m=g.particle_point_m,normal_inward=g.normal_inward))
    for length in [0.,1e-6,4e-6]:
        for height in [.1e-6,.7e-6,3e-6]:
            shape=Capsule([0,0,height],[1,0,2],.5e-6,length);g=triangle_gap(shape,tri)
            exact=height-shape.radius_m-length/np.sqrt(5)
            capsules.append(dict(shape_mode=shape.mode,center_m=shape.center_m,axis=shape.axis_world,R_m=shape.radius_m,L_m=length,
                                 gap_m=g.gap_m,analytic_gap_m=exact,error_m=abs(g.gap_m-exact)))
    return records,ellipsoids,capsules
