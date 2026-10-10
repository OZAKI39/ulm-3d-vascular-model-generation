"""Single predeclared smooth morphology on the native wall; no mechanics search."""
import json, shutil, time
import numpy as np
import pyvista as pv
from .config import ROOT, dump
from .y_bifurcation_sources import YSources, digest, OUT as OLD
from .y_mural_clot_geometry import CONFIG as BASE, build
from .y_surface_coordinates import validate_and_map

TASK='y_bifurcation_irregular_clot_and_initial_fracture_v1'
OUT=ROOT/'verification'/TASK
RESULTS=ROOT/'results'/TASK
EXE=ROOT/'tools/lammps-contact-onset-v1-build/lmp'
CONFIG={**BASE,'geometry_id':'CURVED_MURAL_IRREGULAR_V2','a_m':.00068,'b_m':.00062,'H_m':.00065,
        'geometry_seed':20261008,'randomness_used':False,'irregularity_fraction':.10,
        'morphology_version':'SMOOTH_WARP_THREE_BROAD_LOBES_V1',
        'dt_s':.5e-6,'duration_s':.0011605,'explicit_steps':4660,'dump_every':20,
        'numerical_drag_rate_s_inv':500.,'mpi_ranks':2,'material_id':'LPS_ELASTIC_MAPPING_FROM_LIU2024'}

def raw_profile(x,c):
    p=np.atleast_2d(x);u=(p[:,0]-c['x0_m'])/c['a_m'];v=c['Rv_m']*np.arctan2(p[:,1],p[:,2])/c['b_m']
    U=u-.08*u*u+.04
    V=(v-.07*(1-u*u)-.05*u)/(1+.06*np.sin(np.pi*u))
    q2=U*U+V*V
    # Physical Gaussian sigma >= 0.31 mm: no sub-2-delta features.
    bumps=sum(np.exp(-.5*((u-a)**2+(v-b)**2)/.50**2) for a,b in [(-.40,-.32),(-.18,.40),(.42,.02)])/3
    f=.45*(-np.tanh(1.5*u))+.55*bumps
    return q2,c['H_m']*np.maximum(1-q2,0)**1.5*(1+c['irregularity_fraction']*f)

def profile(x,c):
    q,h=raw_profile(x,c);return q,h/c['height_normalization']

DEFINITION='''# 唯一冻结形态：CURVED_MURAL_IRREGULAR_V2

这是工程合成附壁形态，不是患者 CT/MRI 分割。只在原生 WALL 面片内沿径向向管腔内构造高度。
u=(x+6mm)/a, v=Rv atan2(y,z)/b；a=0.68mm, b=0.62mm, H=0.65mm。
U=u-0.08u²+0.04；V=[v-0.07(1-u²)-0.05u]/[1+0.06 sin(pi u)]；q²=U²+V²。
足迹 q²<1，是平滑弯曲、周向不对称、下游略延长的非椭圆。
B=sum exp(-[(u-u_j)²+(v-v_j)²]/(2*0.50²))/3，中心依次(-.40,-.32),(-.18,.40),(.42,.02)。
f=.45[-tanh(1.5u)]+.55B；|f|<=1；h_raw=H max(1-q²,0)^(3/2)*(1+0.10 f)。
上游高度偏厚，三项宽高斯只形成缓和肩部，不保证三个孤立局部极大值。最小 sigma=0.31mm >2delta=0.301mm。
height_normalization 为预先在固定2001×2001参数格上求得的 max(h_raw/H)，使名义最高高度为0.65mm；不是根据破裂结果调形。
全部高度项随 envelope 在边界归零；无随机噪声。seed仅为任务身份，没有抽样。
沿用真实面片线性细分、全局50um cubic格、27支持点保守管腔掩码与迭代去除六邻居度<2细须的通用规则。
这一个固定几何在所有案例中相同。坐标均为原Y SI全局坐标；不弯曲旧粒子或复用旧连接。

载荷：固定空间 Gaussian中心通过预定参数点u=-.22,v=.04对应实际暴露面最近点确定；sigma=0.28mm。
Opposing 两个中心由u=±.38,v=.08对应实际面确定，sigma=0.28mm，正x/负x两项之差。
每个初始暴露节点使用新V2面积权重与实际面法向；轴向投影 e_x-n_i(n_i·e_x) 不单位化。
每步用当前位置采样空间Gaussian包络，中心不跟随粒子；初始表面切向基冻结在原坐标系，未更新变形表面法向。
因此严格算子为 F_i=A_i P_i^V2 e_x g(x_i(t)) scale pulse(t)：空间场包络为Eulerian，切向投影和面积为参考表面离散算子。
这不是沿变形面重新求解的流体剪切，不对新增断裂面加载。该限制单列在结果中。
0/100/200/100Pa预算为参考态 sum|F|/sum A；所有脉冲1.1605ms、dt0.5us、观察2.33ms。
使用历史exact-membership lineage；persistent沿用至少3个连续保存时刻(本轮观测跨度>=20us)，不是临床持续稳定。
''' 

def freeze():
    assert not (OUT/'GEOMETRY_V2_CONFIG.json').exists()
    c=dict(CONFIG);u,v=np.meshgrid(np.linspace(-1.2,1.2,2001),np.linspace(-1.2,1.2,2001),indexing='ij')
    q=np.c_[c['x0_m']+c['a_m']*u.ravel(),c['Rv_m']*np.sin(c['b_m']*v.ravel()/c['Rv_m']),c['Rv_m']*np.cos(c['b_m']*v.ravel()/c['Rv_m'])]
    c['height_normalization']=float(raw_profile(q,c)[1].max()/c['H_m'])
    dump(OUT/'GEOMETRY_V2_CONFIG.json',c);(OUT/'MORPHOLOGY_DEFINITION.md').write_text(DEFINITION)
    dump(OUT/'MORPHOLOGY_FREEZE.json',dict(frozen_before_geometry_and_mechanics=True,configuration_sha256=digest(OUT/'GEOMETRY_V2_CONFIG.json'),definition_sha256=digest(OUT/'MORPHOLOGY_DEFINITION.md'),time_epoch_s=time.time()))

def geometry():
    assert not (OUT/'GEOMETRY_SUMMARY.json').exists()
    c=json.loads((OUT/'GEOMETRY_V2_CONFIG.json').read_text());s=YSources()
    assert s.mesh.n_points==26712 and s.mesh.n_cells==146230
    g,c,solid,top,attach=build(s,c,out=OUT,results=RESULTS,profile_fn=profile)
    report,clearance=validate_and_map(s,g,c,solid,top,attach,out=OUT)
    # Actual exposed-surface tangents, distinct from the wall tangent used in V1.
    ids=np.flatnonzero(g['surface']);norm=g['normals'][ids];t=np.tile([1.,0.,0.],(len(ids),1));t-=norm[:,0,None]*norm
    assert np.max(np.abs(np.sum(t*norm,axis=1)))<1e-12
    z=dict(np.load(OUT/'SURFACE_NODES_AND_NORMALS.npz'));z.update(axial_tangent_xyz=t/np.linalg.norm(t,axis=1)[:,None],circumferential_tangent_xyz=np.cross(norm,t/np.linalg.norm(t,axis=1)[:,None]),axial_projected_xyz=t,semantics='INITIAL_EXPOSED_SURFACE_QUADRATURE_AND_FROZEN_REFERENCE_TANGENT; CURRENT_POSITION_LOAD_ENVELOPE')
    np.savez_compressed(OUT/'Y_CLOT_V2_SURFACE_NODES_AND_NORMALS.npz',**z)
    np.savez_compressed(OUT/'Y_CLOT_V2_GEOMETRY.npz',**{k:v for k,v in g.items() if k!='info'})
    np.savez_compressed(OUT/'Y_CLOT_V2_BONDS_AND_ANCHORS.npz',bonds=g['bonds'],bond_types=g['bond_types'],anchor_positions=g['X'][g['anchors']],**{k:v for k,v in np.load(OUT/'CURVED_WALL_ANCHORS.npz').items() if k not in ['X']})
    old=json.loads((OLD/'GEOMETRY_SUMMARY.json').read_text());new=g['info'];new['exposed_nodes']=len(ids)
    for key in list(new):
        if key.startswith('old_hemisphere_'):new.pop(key)
    dump(OUT/'GEOMETRY_SUMMARY.json',new)
    diff=new['clot_volume_m3']/new['continuous_mesh_volume_m3']-1
    dump(OUT/'V1_V2_GEOMETRY_COMPARISON.json',dict(V1=old,V2=new,relative_continuous_volume_change=new['continuous_mesh_volume_m3']/old['continuous_mesh_volume_m3']-1,discrete_relative_continuous_error=diff))
    dump(OUT/'GEOMETRY_WALL_CLEARANCE_AND_VOLUME.json',dict(**report,clearance=clearance,discrete_relative_continuous_error=diff,nominal_a_b_H_m=[c[k] for k in ['a_m','b_m','H_m']],shape_frozen=True,surface_attachment_triangles_disjoint=True,dual_topological_node_tags=int((g['surface']&g['base']).sum())))
    # Exact straight WALL chords must also not cross the original lumen.
    from .y_mural_clot_geometry import locator,hit
    loc=locator(s.exterior);bad=sum(hit(loc,g['X'][a],g['X'][b])[0] for a,b in g['bonds'][g['bond_types']==3]);assert bad==0
    dump(OUT/'EXACT_WALL_BOND_CHECK.json',dict(status='PASS',intersecting_wall_chords=bad))
    print('V2_GEOMETRY_GATE_PASS',new,flush=True)

def load_geometry():
    g=dict(np.load(OUT/'Y_CLOT_V2_GEOMETRY.npz'));g['clot_count']=int(g['clot_count']);g['info']=json.loads((OUT/'GEOMETRY_SUMMARY.json').read_text());return g

if __name__=='__main__':
    import sys
    {'freeze':freeze,'geometry':geometry}[sys.argv[1]]()
