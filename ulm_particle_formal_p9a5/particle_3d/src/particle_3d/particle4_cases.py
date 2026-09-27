"""P4 validation-only fixtures, locked population samples and evidence helpers."""
from functools import lru_cache
from itertools import combinations_with_replacement
from dataclasses import asdict
import numpy as np
from .particle_shapes import Sphere,Ellipsoid,Capsule,unit
from .rbc_orientation import rotation_matrix,quaternion_from_short_axis
from .particle3_cases import population_and_geometries,write_json,write_rows,json_safe
from .rbc_capillary_surrogate import area_feasible_interval,capsule_length
from .pair_geometry import pair_gap
from .kinematic_contact import ContactConstraint,project_contacts

P4_BRANCH='dev/particle-4-particle-contact-20260920'
P3_COMMIT='95e54fe474349c35aaad2e2754aff0bad7108c42'
P3_ACCEPTANCE_COMMIT='d37ea16'


@lru_cache(maxsize=1)
def validation_geometries():return population_and_geometries()[1]


def shapes_for_geometry(index=2):
    g=validation_geometries()[index]
    radius_low,radius_high,_=area_feasible_interval(g);radius=(radius_low+radius_high)/2
    return [Sphere(np.zeros(3),.588015722765549e-6),
            Ellipsoid(np.zeros(3),[g.a_m,g.b_m,g.c_m],rotation_matrix(quaternion_from_short_axis([1,2,3]))),
            Capsule(np.zeros(3),unit([1,-1,2]),radius,capsule_length(g.volume_m3,radius))]


def place_pair(a,b,normal,gap=0.):
    n=unit(normal)
    return a,b.moved(a.support(n)-(b.support(-n)-b.center_m)+gap*n)


def shape_record(s):
    return dict(type=type(s).__name__,**asdict(s))


def geometry_cases():
    base=shapes_for_geometry();rows=[]
    for i,j in combinations_with_replacement(range(3),2):
        for direction in [[1,0,0],[1,.3,.7],[-.2,1,.4]]:
            for target in [1e-6,0.,-.02e-6]:
                a,b=place_pair(base[i],base[j],direction,target);g=pair_gap(a,b,17,42);reverse=pair_gap(b,a,42,17)
                rows.append(dict(pair_type=f'{type(a).__name__}-{type(b).__name__}',shape_i=shape_record(a),shape_j=shape_record(b),
                    target_gap_m=target,gap=g.to_dict(),gap_error_m=abs(g.gap_m-target),reciprocity_gap_error_m=abs(g.gap_m-reverse.gap_m),
                    reciprocity_normal_error=np.linalg.norm(g.normal_j_to_i+reverse.normal_j_to_i),
                    reciprocity_point_error_m=max(np.linalg.norm(g.point_i_m-reverse.point_j_m),np.linalg.norm(g.point_j_m-reverse.point_i_m))))
    return rows


def projection_case(shapes,velocities,omegas=None,wall_constraints=()):
    omegas={i:np.zeros(3) for i in shapes} if omegas is None else omegas
    ids=sorted(shapes);gaps=[pair_gap(shapes[i],shapes[j],i,j) for index,i in enumerate(ids) for j in ids[index+1:]]
    constraints=[ContactConstraint.pair(g) for g in gaps if g.state=='TOUCHING']+list(wall_constraints)
    result=project_contacts(shapes,velocities,omegas,constraints)
    return dict(shapes={i:shape_record(s) for i,s in shapes.items()},free_velocities=velocities,free_omegas=omegas,
        gaps=[g.to_dict() for g in gaps],velocities=result.velocities,omegas=result.omegas,
        translation_corrections=result.translation_corrections,angular_corrections=result.angular_corrections,audit=result.record)


def mixed_scene():
    g0,g1=validation_geometries()[0],validation_geometries()[3]
    base=[Sphere(np.zeros(3),.588015722765549e-6),Sphere(np.zeros(3),.7e-6),
          Ellipsoid(np.zeros(3),[g0.a_m,g0.b_m,g0.c_m],np.eye(3)),
          Ellipsoid(np.zeros(3),[g1.a_m,g1.b_m,g1.c_m],np.eye(3))]
    shapes={101:base[0]};ids=[101,203,307,409]
    for k in range(1,4):
        previous=shapes[ids[k-1]]
        shapes[ids[k]]=base[k].moved(previous.support([1,0,0])-base[k].support([-1,0,0])+np.array([(.4+.2*k)*1e-6,0,0]))
    velocities={i:np.array([speed,0,0])*1e-6 for i,speed in zip(ids,[3,1,-1,-3])}
    return shapes,lambda i,s,t:(velocities[i],np.zeros(3))
