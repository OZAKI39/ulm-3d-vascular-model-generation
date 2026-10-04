"""Read-only source/fixture audit before defining any refined physical problem.

The conservative comparison uses the SAME 27 coordinates. It is a separately
labelled proposal, never a replacement or a refinement of the legacy fixture.
"""
from pathlib import Path
import json
import platform
import shutil
import subprocess
import hashlib
import numpy as np
import scipy
import numba
import psutil
from pd_clot.pfpd_candidate.adr_adapter_at1 import fixture, pair_weights
from pd_clot.pfpd_candidate.coupling_stability_audit import physics_config
from pd_clot.pfpd_candidate.at1_evolution import critical_driver
from pd_clot.pfpd_candidate.kernel import Family
from pd_clot.pfpd_candidate.baqp import evaluate

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'verification/pfpd_at1_spatial_discretization_audit_v1'


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def geometry_facts():
    config = physics_config()
    cloud, family, material, top = fixture(config)
    lo, hi = family.X.min(axis=0), family.X.max(axis=0)
    h = cloud.spacing
    origin = np.asarray(config['geometry']['origin_m'])
    cells = np.asarray(config['geometry']['cells'])
    box_hi = origin + cells * h
    counts = np.bincount(family.pairs.ravel(), minlength=len(family.X))
    original = evaluate(family, family.X, material)
    rank = np.linalg.matrix_rank(original.K, tol=None)
    surface = np.any(np.isclose(family.X, lo, rtol=0, atol=1e-15) |
                     np.isclose(family.X, hi, rtol=0, atol=1e-15), axis=1)
    sphere_inside_cells = np.minimum(family.X-origin, box_hi-family.X).min(axis=1) >= family.horizon
    faces_top = cloud.face_normal[:, 2] == 1
    yc = critical_driver(config['phase']['Gc_demo_J_per_m2'], family.horizon)
    facts = dict(
        source_geometry='X = origin + (ijk + 0.5) * h; V_i = h**3',
        source_boundary='fixture overrides fixed mask with minimum/maximum point z planes',
        coordinate_min_m=lo.tolist(), coordinate_max_m=hi.tolist(),
        coordinate_extent_m=(hi-lo).tolist(), spacing_m=h,
        physical_horizon_m=family.horizon, horizon_over_spacing=family.horizon/h,
        generated_cell_domain_min_m=origin.tolist(), generated_cell_domain_max_m=box_hi.tolist(),
        generated_cell_domain_extent_m=(box_hi-origin).tolist(),
        particle_volume_m3=float(cloud.volume[0]), represented_volume_m3=float(cloud.volume.sum()),
        generated_cell_domain_volume_m3=float(np.prod(box_hi-origin)),
        coordinate_hull_volume_m3=float(np.prod(hi-lo)),
        represented_over_hull_volume=float(cloud.volume.sum()/np.prod(hi-lo)),
        generated_top_face_area_m2=float(cloud.face_area[faces_top].sum()),
        coordinate_hull_cross_section_m2=float(np.prod((hi-lo)[:2])),
        bottom_prescribed_point_plane_m=float(lo[2]), top_prescribed_point_plane_m=float(hi[2]),
        physical_clamp_center_separation_m=float(hi[2]-lo[2]),
        voxel_interpretation_clamp_thickness_m=h,
        voxel_interpretation_free_material_gap_m=float(hi[2]-lo[2]-h),
        exact_plane_interpretation_free_length_m=float(hi[2]-lo[2]),
        bottom_points=int((cloud.fixed & ~top).sum()), top_points=int(top.sum()),
        free_points=int((~cloud.fixed).sum()), particles=len(family.X), bonds=len(family.pairs),
        family_counts=counts.tolist(), family_min=int(counts.min()), family_max=int(counts.max()),
        family_mean=float(counts.mean()), family_median=float(np.median(counts)),
        surface_family_mean=float(counts[surface].mean()),
        geometric_non_surface_family_mean=float(counts[~surface].mean()),
        full_untruncated_horizon_points=int(sphere_inside_cells.sum()),
        boundary_truncated_horizon_points=int((~sphere_inside_cells).sum()),
        support_rank_counts={str(r):int((rank == r).sum()) for r in np.unique(rank)},
        Gc_demo_J_m2=config['phase']['Gc_demo_J_per_m2'], Yc_Pa=yc,
        intended_refinement_physical_domain='NOT_UNIQUELY_DECLARED_BY_TENSILE_FIXTURE',
        intended_refinement_clamp_measure='POINT_PLANES_VS_REPRESENTED_VOXEL_SLABS_UNRESOLVED',
    )
    return facts, cloud, family, material, top


def conservative_hull_volumes(X):
    """Tensor trapezoidal weights for an explicitly ASSUMED coordinate-hull box."""
    volume = np.ones(len(X))
    for axis in range(3):
        nodes, indices = np.unique(X[:, axis], return_inverse=True)
        lengths = np.diff(nodes)
        weights = np.r_[lengths[0]/2, (lengths[:-1]+lengths[1:])/2, lengths[-1]/2]
        volume *= weights[indices]
    return volume


def main():
    facts, cloud, family, material, top = geometry_facts()
    write(OUT/'preaudit/GEOMETRY_FACTS.json', facts)
    np.savez_compressed(OUT/'preaudit/LEGACY_FIXTURE.npz', X=family.X, volume=family.volume,
                        pairs=family.pairs, fixed=cloud.fixed, top=top)
    V = conservative_hull_volumes(family.X)
    alternative = Family.from_arrays(family.X, V, family.pairs, family.horizon, family.omega)
    strain = .01
    x = family.X.copy(); x[:, 2] += strain * (family.X[:, 2]-family.X[:, 2].min())
    rows = []
    for name, f in [('LEGACY_FULL_VOXEL_VOLUME', family),
                    ('PROPOSED_HULL_DOMAIN_BOUNDARY_QUADRATURE_NOT_ADOPTED', alternative)]:
        state = evaluate(f, x, material)
        rows.append(dict(case=name, coordinates_unchanged=bool(np.array_equal(f.X, family.X)),
                         volume_m3=float(f.volume.sum()), elastic_energy_J=state.energy,
                         energy_per_represented_volume_Pa=state.energy/f.volume.sum(),
                         energy_per_coordinate_hull_volume_Pa=state.energy/facts['coordinate_hull_volume_m3'],
                         nominal_reaction_N=float(-state.force[top, 2].sum()),
                         free_force_L2_N=float(np.linalg.norm(state.force[~cloud.fixed])),
                         maximum_Fbar_error=float(np.max(abs(state.Fbar-np.diag([1,1,1+strain])))),
                         maximum_Ftilde_error=float(np.max(abs(state.Ftilde-np.diag([1,1,1+strain])))),
                         minimum_shape_rank=int(np.linalg.matrix_rank(state.K).min()),
                         interpretation='Same prescribed affine coordinates; not a new solved/stable branch'))
    write(OUT/'comparisons/ISOLATED_COARSE_QUADRATURE_COMPARISON.json', rows)
    write(OUT/'preaudit/QUADRATURE_CONSISTENCY_GATE.json', dict(
        status='CURRENT_RULE_NOT_PHYSICALLY_CONSISTENT',
        status_scope='Only the proposed fixed 0.25 mm point-span refinement with unchanged V_i=h**3 rule',
        legacy_cell_domain_volume_consistent=True,
        legacy_point_hull_domain_volume_consistent=False,
        physical_problem_selected=False, refinement_generation_authorized_by_gate=False,
        action='NEED_GEOMETRY_DECISION',
        reason='The generator represents a 0.375 mm cell box; tensile clamps are point planes 0.25 mm apart. The continuum domain and clamp measure to preserve jointly are not declared.',
        diagnostic_refinement_volume_formula='For fixed point span L and n points: sum V = L**3 * (n/(n-1))**3; not constant.',
        proposed_corrected_coarse_scope='A 0.25 mm hull-domain volume comparison only, NOT adopted; no refined lattice generated.'))
    mem = psutil.virtual_memory()
    hw = dict(cpu=platform.processor(), platform=platform.platform(),
              physical_cores=psutil.cpu_count(logical=False), logical_cores=psutil.cpu_count(),
              RAM_total_bytes=mem.total, RAM_available_bytes=mem.available,
              disk_free_bytes=shutil.disk_usage(ROOT).free, Python=platform.python_version(),
              NumPy=np.__version__, SciPy=scipy.__version__, Numba=numba.__version__,
              reserve_total_RAM_fraction=.2, dense_available_RAM_max_fraction=.4,
              Julia='Not launched: geometry hard gate precedes external mechanics comparisons')
    write(OUT/'provenance/HARDWARE.json', hw)
    sources=['pd_clot/geometry.py','pd_clot/pfpd_candidate/adr_adapter_at1.py',
             'pd_clot/pfpd_candidate/kernel.py','pd_clot/pfpd_candidate/baqp.py',
             'verification/pfpd_at1_evolution_v1/provenance/PREDECLARED_CONFIG.json']
    write(OUT/'preaudit/SOURCE_IDENTITIES.json', {s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources})
    report = f'''# Geometry and quadrature hard-gate audit

状态：**NEED_GEOMETRY_DECISION**。尚未生成任何 refined lattice。

## 源码确定的离散事实

生成器 `pd_clot/geometry.py` 使用 `X=origin+(ijk+0.5)h`、`V_i=h³`，同时把外表面放在最外层点外侧半个 h。
张拉适配器 `pd_clot/pfpd_candidate/adr_adapter_at1.py` 则直接约束最小/最大点坐标平面；`grip_center_separation_m` 元数据是0.25 mm，并没有定义整个材料域边长。

|量|实际值|
|---|---:|
|点坐标范围（三轴相同）|−0.125 至 +0.125 mm|
|点跨度／夹持中心距离|0.25 mm|
|生成器体素域（三轴相同）|−0.1875 至 +0.1875 mm|
|生成器体素域边长|0.375 mm|
|点间距 h|0.125 mm|
|固定物理视界 δ|0.25125 mm|
|δ/h|2.01|
|每点体积|0.001953125 mm³|
|总表示体积|0.052734375 mm³|
|点包围盒体积|0.015625 mm³|
|体积比|3.375|
|生成器完整顶面面积|0.140625 mm²|
|点包围盒截面积|0.0625 mm²|
|点数／键数|27／185|
|顶／底约束点数|9／9|
|自由点数／自由位移分量|9／27|
|family min／median／mean／max|{facts['family_min']}／{facts['family_median']}／{facts['family_mean']:.12g}／{facts['family_max']}|
|满球形邻域／截断邻域点数|{facts['full_untruncated_horizon_points']}／{facts['boundary_truncated_horizon_points']}|
|shape rank|全部3|

0.375 mm并不是凭空推测的尺寸，它来自生成器的体素边界和体积。0.25 mm也不是错误的夹持距离，它来自实际施加位移的点平面。**未确定的是细化研究要固定哪一种材料域和夹持区域的连续解释。**

## 为什么直接加密会改变物理问题

如果把0.25 mm点跨度固定，用n点和h=L/(n−1)，但继续让每点代表完整h³，则总量是 `L³[n/(n−1)]³`，随n变化。
如果固定生成器的0.375 mm体素域，则n点的h=L/n；最外点夹持间距变为L−h，原先“夹持最外一层”的物理位置和厚度都会变化。
因此两种直接复用方式都不能同时保持提示词要求的物理域、体积和夹持定义。

## 夹持条件的未决部分

当前代码准确规定z=±0.125 mm处的9个点固定，但未声明其连续含义是：

- 0.25 mm材料域的上下边界面；或
- 0.375 mm体素材料域中的上下整层刚性夹持材料。按完整体素解释，两夹持层各厚0.125 mm，内侧间自由材料厚0.125 mm；代码没有给出跨分辨率延拓。

这会改变边界面积、自由长度及允许位移，不应由数值程序静默选择。

## 已完成的隔离比较

保持27点坐标、185键、δ、材料、Gc完全不变，另计算了**尚未采用**的0.25 mm点包围盒张量梯形体积权重。面／棱／角点分别按边界对偶体积积分，总量为0.015625 mm³。
原始规则和提议规则均对同一仿射位形计算能量、力及F，结果见 `../comparisons/ISOLATED_COARSE_QUADRATURE_COMPARISON.json`。
该比较既不是新的已验证粗网格，也不是旧分支的逐位延续，未创建新网格和新稳定性结论。

## 需要明确的几何选择

建议：把本轮基准定义为**0.25 mm边界节点立方体**，上下z=±0.125 mm为精确夹持面，使用守恒边界体积权重。旧27点完整体素结果独立保留为legacy；修正后的coarse单独命名并重新通过mechanics／threshold检查，然后才开展5³、7³。

另一种定义：坚持**0.375 mm完整体素域**，显式规定原0.125 mm厚的上下刚性夹持区域及其运动，再另行确定能保持该定义的细化点云。这需要改变预期的点布局，不能沿用0.25 mm节点边界方案。

停止依据是本轮用户提示词第3节的hard gate，以及第4／6节禁止静默改变物理边界和积分规则。该停止不来自SKILL.md或额外审批规则。
'''
    (OUT/'preaudit/GEOMETRY_AND_QUADRATURE_AUDIT.md').write_text(report)
    print(json.dumps(dict(facts=facts,coarse_comparison=rows),indent=2))


if __name__ == '__main__':
    main()
