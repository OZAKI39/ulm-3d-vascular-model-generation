# Open inlet 审核：PASS

WALL 只包含真实 WALL。INLET cap 179 个面与 WALL 的 GlobalNodeID 面集合交集为 0；perimeter 顶点属于真实 WALL，保留 rim 碰撞。面积 7.756795804427715e-12 m²，等效直径 3.1426516127 µm。cap 略非平面，最大平面残差 1.48136e-11 m；查询使用实际三角和 owner tetra，不把拟合平面当碰撞面。

永久 tests 已运行并通过：

| 验证 | 结果 |
| --- | --- |
| sphere_straddles_open_cap_accepts | 开放圆管入口允许同一球向上游跨出 |
| same_sphere_with_solid_cap_rejects | 加入实体 cap 的对照拒绝 |
| rim_wall_overlap_rejects | 真实 rim/WALL overlap 拒绝 |
| inlet_cap_not_in_wall_bvh | 实际 cap/WALL 面交集 0 |
| center_on_inlet_field_sample_is_valid | 179 cap 面心均有效 |
| edge_and_vertex_boundary_ownership_stable | 重复及反序查询 owner 稳定 |

上述测试分布于 `test_open_inlet_straddling_sphere_accepts.py`、`test_real_wall_overlap_rejects.py`、`test_inlet_cap_not_solid.py` 和旧 open-cap 测试。真实100k proposals 没有 CENTER_OUTSIDE 拒绝；WALL 80672、lower handoff 107。

birth center 位于真实 cap，只要求 center 有效及球不碰真实 WALL/不越过 lower handoff；不要求整球 inside，也没有 distance_to_inlet_plane≥radius。较大尺寸接近零 acceptance 来源于狭窄真实壁面，不是把 cap 封死。

30 smoke 全部保存线段按向外方向与真实 INLET 相交事后复判：逃逸 0。初始 cap 上向内运动不误记为逃逸。此观察只证明本次保存轨迹与原连续 WALL 证书的结果；不构成任意上游球扫掠几何的证明，也没有新增 upstream-to-cap 动力学。

[几何机器证据](data/open_inlet_geometry.json) · [最终测试日志](logs/final_new_and_legacy_pass.txt) · [smoke 逐轨迹指标](data/smoke30_metrics.json)
