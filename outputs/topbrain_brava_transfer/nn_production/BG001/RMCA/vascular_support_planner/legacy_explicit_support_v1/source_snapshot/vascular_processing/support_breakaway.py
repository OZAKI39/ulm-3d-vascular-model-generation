"""Controlled bottom-only feet and physical coupon geometry, never strength claims."""
from __future__ import annotations
import numpy as np
import trimesh
from .support_tree_planner import strut
from .support_inputs import save,table

def widths(profile):
    w=max(profile['line_width'],profile['nozzle_diameter'])
    return [round(w*n,6) for n in (2,3,4)]

def root_geometry(xy,floor,trunk_diameter,neck_width,profile,root_id='R0',load=1):
    h=profile['layer_height'];foot_h=4*h;neck_h=3*h
    foot_d=max(trunk_diameter+4*profile['line_width'],neck_width*2)
    z0=floor-.05;z1=floor+foot_h;z2=z1+neck_h
    # Upward widening after the neck remains below 35 degrees from vertical.
    z3=z2+(trunk_diameter-neck_width)/2/np.tan(np.radians(35))
    at=lambda z:np.r_[xy,z]
    foot=strut(at(z0),at(z1),foot_d)
    neck=strut(at(z1-.02),at(z2+.02),neck_width)
    flare=strut(at(z2),at(z3),neck_width,trunk_diameter)
    mesh=trimesh.boolean.union([foot,neck,flare],engine='manifold')
    record=dict(root_id=root_id,xyz=at(floor).tolist(),contact_area_mm2=np.pi*foot_d**2/4,
        foot_diameter_mm=foot_d,neck_width_mm=neck_width,neck_height_mm=neck_h,neck_lines=neck_width/profile['line_width'],
        neck_bottom_mm=z1,neck_top_mm=z2,trunk_start_mm=z3,trunk_diameter_mm=trunk_diameter,trunk_load=load,
        estimated_break_location=at((z1+z2)/2).tolist(),root_role='APPROVED_INNER_BOTTOM_BREAKAWAY',break_force_validated=False)
    return mesh,foot,record

def bottom_root_allowed(root,inner,margin=.05):
    lo,hi=np.asarray(inner);p=np.asarray(root['xyz']);r=root['foot_diameter_mm']/2
    return bool(abs(p[2]-lo[2])<1e-6 and np.all(p[:2]-r>=lo[:2]+margin) and np.all(p[:2]+r<=hi[:2]-margin))

def coupon(out,profile):
    out.mkdir(parents=True,exist_ok=True);parts=[];rows=[];features=[];h=profile['layer_height']
    for a,neck in enumerate(widths(profile)):
        for b,tip in enumerate(widths(profile)):
            for c,gap in enumerate([.15,.20,.25]):
                n=a*9+b*3+c;xy=np.array([35+14*(n%9),35+12*(n//9)],float);floor=2.;rod_z=10.;rod_r=.6
                slab=trimesh.creation.box([12,10,2]);slab.apply_translation(np.r_[xy,1.])
                root,foot,r=root_geometry(xy,floor,2.2,neck,profile,f'COUPON_{n:02}')
                tiptop=np.r_[xy,rod_z-rod_r-gap];shoulder=tiptop-[0,0,1.6]
                trunk=strut(np.r_[xy,r['trunk_start_mm']-.02],shoulder+[0,0,.02],2.2,1.2)
                taper=strut(shoulder,tiptop,1.2,tip)
                piers=[strut(np.r_[xy+[x,0],floor-.02],np.r_[xy+[x,0],rod_z],1.6) for x in [-4.,4.]]
                rod=strut(np.r_[xy+[-4.5,0],rod_z],np.r_[xy+[4.5,0],rod_z],2*rod_r)
                parts.append(trimesh.boolean.union([slab,root,trunk,taper,*piers,rod],engine='manifold'))
                rows.append(dict(coupon_id=f'COUPON_{n:02}',row=n//9,column=n%9,center_xy_mm=xy.tolist(),neck_width_mm=neck,tip_width_mm=tip,top_gap_mm=gap,physical_test_completed=False))
                features.extend([dict(feature_id=f'{n}_TIP',kind='tip',xyz=(tiptop-[0,0,h/2]).tolist(),width_mm=tip),
                    dict(feature_id=f'{n}_NECK',kind='neck',xyz=r['estimated_break_location'],width_mm=neck)])
    from .support_custom_qc import canonical_stl
    mesh=canonical_stl(trimesh.util.concatenate(parts));mesh.export(out/'support_breakaway_coupon.stl')
    table(out/'coupon_layout.csv',rows);save(out/'coupon_features.json',features)
    save(out/'coupon_manifest.json',dict(cases=rows,tip_widths_mm=widths(profile),neck_widths_mm=widths(profile),top_gaps_mm=[.15,.20,.25],physical_validation=False))
    (out/'coupon_test_plan.md').write_text('''# 折断根部与血管接触端试验

先打印本试片，再决定是否打印完整血管件。它包含 27 个独立单元：三种颈宽 × 三种尖端宽度 × 三种顶部间隙，排列位置见 coupon_layout.csv。每格均有 2 mm 底板、可折断支撑、小直径水平 ABS 圆杆和两侧永久支柱；不是简单立方体。打印工艺沿用 X1 Carbon 0.2 mm、ABS Vessel，关闭自动支撑。

1. 冷却后先拍照，记录每格颈部、支撑尖端是否完整，圆杆是否下垂或被喷嘴碰撞。
2. 从顶部夹住支撑，记录在颈部折断还是其它部位破坏；如有测力计，同时记录施力方向与峰值力。软件不预测合格折断力。
3. 拆除尖端，记录圆杆是否断裂、表面是否拉伤、间隙是否黏连，以及实际拆除难易。
4. 每种组合至少重复三次；记录打印温度、冷却时间、材料批次和照片。优先选择可稳定打印、在预定颈部折断且圆杆损伤少的组合。

0.15/0.20/0.25 mm 是几何间隙；0.10 mm 层高下实际刀路间隙会量化，必须结合预览记录。不要仅按 STL 间隙判断拆除效果。试片通过不等同于完整血管模型或长短桥已经通过。
''',encoding='utf-8')
    return mesh,rows,features
