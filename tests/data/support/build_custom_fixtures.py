"""Rebuild permanent synthetic explicit-support fixtures; source geometry is unrelated to BG001."""
from pathlib import Path
import json
import numpy as np
import trimesh
from vascular_processing.support_tree_planner import strut
from vascular_processing.support_breakaway import root_geometry

HERE=Path(__file__).resolve().parent
PROFILE=dict(layer_height=.1,first_layer_height=.12,line_width=.22,nozzle_diameter=.2,support_top_gap=.2,support_xy_gap=.45)
def chain(points,r=.6):
    parts=[strut(a,b,2*r) for a,b in zip(points[:-1],points[1:])]
    for p in points[1:-1]:
        ball=trimesh.creation.icosphere(subdivisions=2,radius=r);ball.apply_translation(p);parts.append(ball)
    return trimesh.boolean.union(parts,engine='manifold')
def main():
    t=np.linspace(0,30,41);s=chain(np.c_[t,3*np.sin(t/5),np.full_like(t,10.)]);curve=chain(np.c_[t,np.sin(t/10)*8,np.full_like(t,10.)])
    fork=trimesh.boolean.union([strut([0,0,10],[10,0,10],1.2),strut([9.5,0,10],[20,5,10],1.2),strut([9.5,0,10],[20,-5,10],1.2)],engine='manifold')
    double=trimesh.boolean.union([fork,strut([19,4.5,10],[28,8,12],1.2),strut([19,4.5,10],[28,2,10],1.2)],engine='manifold')
    close=trimesh.util.concatenate([strut([0,0,10],[20,0,10],1),strut([0,1.5,10],[20,1.5,10],1)])
    unstable=trimesh.boolean.union([strut([0,0,0],[0,0,30],.5),strut([0,0,29.9],[12,0,35],.5)],engine='manifold')
    bridge=trimesh.boolean.union([strut([-2,0,0],[-2,0,6],2.2),strut([2,0,0],[2,0,6],2.2),strut([-2,0,5.5],[2,0,5.5],1.1)],engine='manifold')
    root,_,_=root_geometry([0,0],2,2.2,.66,PROFILE)
    torus=trimesh.creation.torus(major_radius=3,minor_radius=.5);torus.apply_translation([0,0,8]);trapped=trimesh.boolean.union([torus,strut([3,0,0],[3,0,8],1)],engine='manifold')
    models={'long_curved_tube':curve,'S_curve':s,'single_bifurcation':fork,'double_bifurcation':double,'thin_vessel':strut([0,0,10],[20,0,10],.7),'close_parallel_vessels':close,'unstable_tree':unstable,'pillar_bridge':bridge,'breakaway_root':root,'trapped_tree':trapped}
    for name,mesh in models.items():
        mesh.export(HERE/(name+'.stl'));(HERE/(name+'.json')).write_text(json.dumps(dict(name=name,synthetic=True,units='mm',volume_mm3=float(mesh.volume),physical_test_completed=False),indent=2)+'\n')
if __name__=='__main__':main()
