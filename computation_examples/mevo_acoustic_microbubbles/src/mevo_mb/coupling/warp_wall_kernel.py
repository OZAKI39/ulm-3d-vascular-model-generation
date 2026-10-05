"""Warp-only module; importing the standard package never requires Warp."""
import warp as wp


@wp.kernel
def query_mesh(mesh:wp.uint64,points:wp.array(dtype=wp.vec3),
               distance:wp.array(dtype=float),closest:wp.array(dtype=wp.vec3),
               normal:wp.array(dtype=wp.vec3),face:wp.array(dtype=int)):
    i=wp.tid()
    point=points[i]
    # A geometric search bound larger than the whole millimetre-scale vessel;
    # this is not a physical force/contact cutoff.
    result=wp.mesh_query_point(mesh,point,1.0)
    face[i]=-1
    if result.result:
        p=wp.mesh_eval_position(mesh,result.face,result.u,result.v)
        delta=point-p
        d=wp.length(delta)
        n=wp.vec3(0.0,0.0,0.0)
        if d>0.0:
            n=delta/d*result.sign
        else:
            n=wp.mesh_eval_face_normal(mesh,result.face)
        distance[i]=d*result.sign
        closest[i]=p
        normal[i]=n
        face[i]=result.face
