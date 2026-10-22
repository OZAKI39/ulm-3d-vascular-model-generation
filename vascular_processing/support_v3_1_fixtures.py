"""Small deterministic five-wall fixtures for orientation/release regression."""
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from .support_v3_geometry import union


def block(size,center):
    return trimesh.creation.box(size,transform=trimesh.transformations.translation_matrix(center))


def inverted_box():
    return union([block([12,12,1],[0,0,5.5]),block([1,12,5],[-5.5,0,2.5]),block([1,12,5],[5.5,0,2.5]),
                  block([10,1,5],[0,-5.5,2.5]),block([10,1,5],[0,5.5,2.5])])


def tilt_fixture():
    mold=union([inverted_box(),block([.4,1,.3],[0,0,1.2])])
    support=union([block([.4,.5,2.2],[1.6,0,1.1]),block([1.7,.5,.2],[.95,0,2.1])])
    return mold,support


def fixture(name):
    box=inverted_box()
    if name=='open_box_inverted':return box
    if name=='open_box_inverted_tilt':
        m=np.eye(4);m[:3,:3]=Rotation.from_euler('y',15,degrees=True).as_matrix();box.apply_transform(m);box.apply_translation([0,0,-box.bounds[0,2]]);return box
    if name=='tube_inside_inverted_box':
        tube=trimesh.creation.cylinder(radius=.4,height=6,sections=32);tube.apply_transform(trimesh.geometry.align_vectors([0,0,1],[1,0,0]));tube.apply_translation([0,0,3]);return union([box,tube])
    if name=='ceiling_failure_box':return box
    if name=='plate_release_pass':return block([.66,3,4],[0,0,2])
    if name=='plate_release_fail':return union([block([.5,.5,4],[4,0,2]),block([8,.5,.4],[0,0,3.8])])
    raise ValueError(name)


def generate(out):
    out.mkdir(parents=True,exist_ok=True)
    names=['open_box_inverted','open_box_inverted_tilt','tube_inside_inverted_box','ceiling_failure_box','plate_release_pass','plate_release_fail']
    for name in names:fixture(name).export(out/(name+'.stl'))
