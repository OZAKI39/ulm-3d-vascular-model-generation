"""Quantized FDM solids: horizontal layer sections, roofs and break sections."""
from __future__ import annotations
import numpy as np
import trimesh
from shapely.geometry import Polygon
from functools import lru_cache

KINDS=('VERTICAL_PILLAR','SELF_SUPPORTING_BRANCH','Y_MERGE','SHORT_BRIDGE','LONGITUDINAL_RAIL','ROOF_PAD','BUTTRESS','BREAKAWAY_FOOT','BREAKAWAY_SEGMENT_NECK')

def quantize(value,unit,minimum=1):return float(max(minimum,np.ceil(value/unit-1e-9))*unit)

def loft(a,b,size_a,size_b=None,angle=0.):
    a=np.asarray(a,float);b=np.asarray(b,float);sa=np.broadcast_to(size_a,(2,));sb=sa if size_b is None else np.broadcast_to(size_b,(2,))
    theta=float(angle);r=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    square=np.array([[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]])
    vertices=np.vstack([np.c_[(square*sa)@r.T+a[:2],np.full(4,a[2])],np.c_[(square*sb)@r.T+b[:2],np.full(4,b[2])]])
    faces=[[0,2,1],[0,3,2],[4,5,6],[4,6,7]]
    for i in range(4):j=(i+1)%4;faces.extend([[i,j,j+4],[i,j+4,i+4]])
    # Eight unique vertices and twelve consistently oriented triangles are
    # constructed explicitly. General-purpose welding is redundant here.
    mesh=trimesh.Trimesh(vertices,faces,process=False)
    mesh._cache['is_convex']=True
    return mesh

def roof_pad(xyz,length,width,thickness,gap,angle=0.):
    xyz=np.array(xyz,float);top=xyz-[0,0,gap]
    mesh=loft(top-[0,0,thickness],top,[length,width],angle=angle)
    return mesh,Polygon(mesh.vertices[mesh.vertices[:,2]>top[2]-1e-7,:2]).convex_hull

def pillar(a,b,width):return loft(a,b,width)

def branch(a,b,width):return loft(a,b,width)

def y_merge(left,right,merge,root,width):
    return trimesh.boolean.union([branch(merge,left,width),branch(merge,right,width),pillar(root,merge,width)],engine='manifold')

def horizontal_member(a,b,width,height):
    a=np.asarray(a);b=np.asarray(b);axis=b[:2]-a[:2];length=np.linalg.norm(axis)
    center=(a+b)/2;angle=np.arctan2(axis[1],axis[0]);return loft(center-[0,0,height/2],center+[0,0,height/2],[length+width,width],angle=angle)

def buttress(bottom,top,bottom_size,top_size,angle=0.):return loft(bottom,top,bottom_size,top_size,angle)

def neck(xyz,width,height,angle=0.):return loft(np.asarray(xyz)-[0,0,height/2],np.asarray(xyz)+[0,0,height/2],width,angle=angle)

@lru_cache(maxsize=256)
def _foot_template(width,h,w,bed):
    foot_width=quantize(width+4*w,w);neck_width=3*w;base=np.zeros(3)
    foot=loft(base if bed else base-[0,0,.03],base+[0,0,4*h],foot_width)
    cut=base+[0,0,5.5*h];neckmesh=neck(cut,neck_width,3*h)
    flare_height=quantize(max(width-neck_width,0)/2/.65,h)
    start=base+[0,0,7*h];end=start+[0,0,flare_height]
    upper=trimesh.boolean.union([neckmesh,buttress(start,end,neck_width,width)],engine='manifold')
    mesh=trimesh.boolean.union([foot,upper],engine='manifold')
    return mesh,foot,cut,end,upper

def breakaway_foot(xy,z,width,profile,domain,identity):
    h=profile['layer_height'];w=profile['line_width'];foot_width=quantize(width+4*w,w);neck_width=3*w
    base=np.r_[xy,z];template,foot_template,cut0,end0,_=_foot_template(width,h,w,domain=='BUILD_PLATE')
    mesh=template.copy();mesh.apply_translation(base);foot=foot_template.copy();foot.apply_translation(base);cut=cut0+base;end=end0+base
    record=dict(root_id=identity,xyz=base.tolist(),domain=domain,foot_width_mm=foot_width,neck_width_mm=neck_width,neck_height_mm=3*h,
        neck_xyz=cut.tolist(),trunk_start=end.tolist(),width_mm=width,approved_box_contact=domain!='BUILD_PLATE',first_layer_area_mm2=foot_width**2 if domain=='BUILD_PLATE' else None)
    return mesh,foot,record

def foot_above_approved_pad(xy,z,width,profile,domain):
    upper=_foot_template(width,profile['layer_height'],profile['line_width'],domain=='BUILD_PLATE')[4].copy()
    upper.apply_translation([*xy,z]);return upper

def primitive_gallery(profile):
    w=profile['line_width'];h=profile['layer_height'];items={}
    items['VERTICAL_PILLAR']=pillar([0,0,0],[0,0,8],8*w)
    items['SELF_SUPPORTING_BRANCH']=branch([0,0,0],[2,0,8],6*w)
    items['Y_MERGE']=y_merge([-2,0,8],[2,0,8],[0,0,3],[0,0,0],6*w)
    items['SHORT_BRIDGE']=trimesh.boolean.union([pillar([-1,0,0],[-1,0,5],6*w),pillar([1,0,0],[1,0,5],6*w),horizontal_member([-1,0,5],[1,0,5],6*w,4*h)],engine='manifold')
    items['LONGITUDINAL_RAIL']=items['SHORT_BRIDGE'].copy()
    items['ROOF_PAD']=roof_pad([0,0,1],8*w,4*w,3*h,0)[0]
    items['BUTTRESS']=buttress([0,0,0],[0,0,3],3*w,8*w)
    items['BREAKAWAY_FOOT']=breakaway_foot([0,0],0,8*w,profile,'BUILD_PLATE','GALLERY')[0]
    items['BREAKAWAY_SEGMENT_NECK']=neck([0,0,.3],3*w,3*h)
    return items
