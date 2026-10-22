"""Reproducible FDM V2 fixtures; all dimensions are mm, no patient data."""
from pathlib import Path
import json,runpy
import numpy as np
import trimesh
from vascular_processing import support_v2_primitives as p
from vascular_processing.support_v2_mesh import support_stl
from vascular_processing.support_tree_planner import strut

HERE=Path(__file__).resolve().parent

def main():
    out=HERE/'v2';out.mkdir(exist_ok=True)
    chain=runpy.run_path(str(HERE/'build_custom_fixtures.py'))['chain']
    x=np.linspace(0,24,33)
    fork=trimesh.boolean.union([strut([0,0,10],[12,0,10],1.2),strut([11.5,0,10],[20,6,10],1.2),strut([11.5,0,10],[20,-6,10],1.2)],engine='manifold')
    straight=strut([0,0,10],[24,0,10],1.2)
    wall=trimesh.creation.box([2,8,8]);wall.apply_translation([12,0,4])
    pillars=[p.pillar([q,0,0],[q,0,8],1.76) for q in [9.7,13.7]]
    rail=p.horizontal_member([9.7,0,8],[13.7,0,8],1.76,.4)
    closed=trimesh.boolean.difference([trimesh.creation.box([10,10,10]),trimesh.creation.box([8,8,8])],engine='manifold')
    models=dict(straight_horizontal_tube=straight,long_curved_tube=chain(np.c_[x,5*np.sin(x/15),np.full_like(x,10.)]),
        S_tube=chain(np.c_[x,3*np.sin(x/4),np.full_like(x,10.)]),Y_bifurcation=fork,
        high_arch=chain(np.c_[x,np.zeros_like(x),4+10*np.sin(np.pi*x/24)]),low_arch=chain(np.c_[x,np.zeros_like(x),4+2*np.sin(np.pi*x/24)]),
        external_port=strut([-8,0,10],[8,0,10],3.),tube_over_wall=trimesh.util.concatenate([straight,wall]),
        rail_supported_tube=trimesh.boolean.union([*pillars,rail],engine='manifold'),
        bridge_supported_tube=trimesh.boolean.union([*pillars,rail],engine='manifold'),trapped_scaffold=closed)
    for name,m in models.items():
        if name=='rail_supported_tube':
            pads=[p.roof_pad([q,0,9.2],1.32,.88,.3,0)[0] for q in [10.,11.7,13.4]]
            posts=[p.pillar([q,0,8],[q,0,8.9],.88) for q in [10.,11.7,13.4]]
            m=trimesh.boolean.union([m,*pads,*posts],engine='manifold')
        if name in ['rail_supported_tube','bridge_supported_tube']:
            m=trimesh.util.concatenate([m,straight]);straight.export(out/(name+'_vascular.stl'))
        m=support_stl(m);m.export(out/(name+'.stl'))
        (out/(name+'.json')).write_text(json.dumps(dict(name=name,synthetic=True,units='mm',volume_mm3=float(m.volume),physical_test_completed=False),indent=2)+'\n')
    (out/'README.md').write_text('# V2 合成几何\n\n运行 `runpy.run_path("tests/data/support/build_v2_fixtures.py", run_name="__main__")` 重建。血管曲线、端口、盒壁、桥轨和封闭腔分别用于欠支撑与拆除反例；不是通过制造验证的模型。*_vascular.stl 是两个支撑场景中独立保存的血管。\n')

if __name__=='__main__':main()
