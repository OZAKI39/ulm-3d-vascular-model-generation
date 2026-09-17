"""
候选分支约束检查。

DCCO 不是“随机连线”算法。每次新增血管都要满足一组几何和生理约束：
- 血管段不能太短，否则会产生小毛刺；
- 分叉角不能太尖或太反向，否则不符合血管形态；
- 线段不能穿过禁区或与已有血管相交；
- 半径更新后仍要满足长径比和 Murray 对称性。

这个模块只负责判断“能不能接受”，不负责真正修改树结构。
"""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np

from .geometry import _angle_deg, _distance
from ..core.tree_index import segment_distances_3d


class ConstraintMixin:
    """把所有接受/拒绝候选分支的规则集中在一起。"""

    def _passes_local_constraints(self, vj_id: int, xb: np.ndarray,
                                   xnew: np.ndarray) -> bool:
        """
        检查普通 split 模式下的局部约束。

        `vj` 是候选父血管，`xb` 是候选分叉点，`xnew` 是新终端点。
        如果返回 False，说明这个候选不应该进入更昂贵的试插入流程。
        """
        # vj 是候选父血管。
        vj = self.vessels[vj_id]
        # non_branching 不允许分叉；distal 模式不走 split 约束。
        if vj.branching_mode in {"non_branching", "distal"}:
            return False
        # 分裂后近端段 vp：vj.x_p -> xb。
        vp_xp, vp_xd = vj.x_p, xb
        # 分裂后远端延续段 vs：xb -> vj.x_d。
        vs_xp, vs_xd = xb, vj.x_d
        # 新终端分支 vnew：xb -> xnew。
        vnew_xp, vnew_xd = xb, xnew

        # 计算三条新段长度。
        lp = _distance(vp_xd, vp_xp)
        ls = _distance(vs_xd, vs_xp)
        lnew = _distance(vnew_xd, vnew_xp)
        # 任一段几乎为 0 都会导致退化几何。
        if lp < 1e-6 or ls < 1e-6 or lnew < 1e-6:
            return False

        # 防毛刺：如果 xb 太靠近父血管任一端点，分裂后会产生非常短的小段。
        # 这些小段在可视化中会像“血管头上的小刺”，在仿真中也会造成不必要
        # 的高阻力段，因此直接拒绝。
        # vj_length 用于按父段比例限制 xb 不能太靠近任一端。
        vj_length = _distance(vj.x_d, vj.x_p)
        min_parent_split_segment = max(
            self.cfg.min_segment_length_um,
            self.cfg.min_split_fraction * vj_length,
        )
        # 父段分裂后的两段都必须足够长。
        if lp < min_parent_split_segment or ls < min_parent_split_segment:
            return False
        # 新终端段也不能短于最小段长。
        if lnew < self.cfg.min_segment_length_um:
            return False

        # 防长尖刺：新终端分支过长通常说明这个终端点离树太远，连接代价高
        # 且形态不自然。
        if self.cfg.max_terminal_segment_um > 0 and lnew > self.cfg.max_terminal_segment_um:
            return False

        # 分叉角约束：vp_dir 是血流到达 xb 的方向，vs_dir/vnew_dir 是两个
        # 子分支离开 xb 的方向。我们拒绝反向折返、过小夹角和过大开口角。
        # vp_dir 表示血流到达 xb 的方向。
        vp_dir = vp_xd - vp_xp
        # vs_dir 表示旧远端延续段离开 xb 的方向。
        vs_dir = vs_xd - vs_xp
        # vnew_dir 表示新分支离开 xb 的方向。
        vnew_dir = vnew_xd - vnew_xp
        # bend_vs/bend_vnew 表示相对父方向的弯折角。
        bend_vs = _angle_deg(vp_dir, vs_dir)
        bend_vnew = _angle_deg(vp_dir, vnew_dir)
        # 子段相对父方向弯折不能过大。
        if (bend_vs > self.cfg.max_daughter_bend_deg
                or bend_vnew > self.cfg.max_daughter_bend_deg):
            return False
        # 新分支也不能几乎沿父方向贴着走。
        if bend_vnew < self.cfg.min_new_branch_angle_deg:
            return False
        # opening 是两个子段之间的开口角。
        opening = _angle_deg(vs_dir, vnew_dir)
        # 开口角太小或太大都拒绝。
        if opening < self.cfg.min_daughter_angle_deg or opening > self.cfg.max_opening_angle_deg:
            return False

        # 分裂 vj 会改变旧孩子看到的“父方向”。旧孩子虽然不是新生成的，
        # 但它们所在 junction 的角度也可能因此变坏，所以必须重新检查。
        for c_id in vj.children:
            child = self.vessels[c_id]
            child_dir = child.x_d - child.x_p
            if _angle_deg(vs_dir, child_dir) > self.cfg.max_daughter_bend_deg:
                return False

        # 分裂 vj 也会改变它在自己父节点看来时的方向。
        # 因此要重新检查上游兄弟分支开口角，避免旧 junction 被间接破坏。
        if vj.parent_id >= 0:
            parent = self.vessels[vj.parent_id]
            for sibling_id in parent.children:
                if sibling_id == vj_id:
                    continue
                sibling = self.vessels[sibling_id]
                sibling_dir = sibling.x_d - sibling.x_p
                sibling_opening = _angle_deg(vp_dir, sibling_dir)
                # vj 新方向和兄弟分支之间也必须满足 opening angle 约束。
                if (sibling_opening < self.cfg.min_daughter_angle_deg
                        or sibling_opening > self.cfg.max_opening_angle_deg):
                    return False
            parent_dir = parent.x_d - parent.x_p
            # vj 新近端方向相对上游父段不能弯折过大。
            if _angle_deg(parent_dir, vp_dir) > self.cfg.max_daughter_bend_deg:
                return False

        # 区域约束：distribution 血管必须留在 perfusion domain；transport 和
        # perforator 可以使用 carriage space，但仍不能进入 avascular region。
        # 分叉点和新终端点必须位于各自允许区域。
        if not (self._bifurcation_domain_ok(vj, xb) and self._in_perfusion_domain(xnew)):
            return False
        # 分裂后的近端段必须符合父血管角色区域约束。
        if not self._segment_domain_ok(vp_xp, vp_xd, vj.role):
            return False
        # 分裂后的远端延续段也必须符合父血管角色区域约束。
        if not self._segment_domain_ok(vs_xp, vs_xd, vj.role):
            return False
        # 新分支按默认血管角色检查区域。
        if not self._segment_domain_ok(vnew_xp, vnew_xd, self.cfg.default_vessel_role):
            return False

        # 碰撞约束：新线段不能穿过已有血管。需要排除共享端点的相邻血管，
        # 否则合法连接会被误判为“相交”。
        # vp 与原 vj 共享几何位置，检查相交时必须排除 vj 自身。
        exclude_vp = {vj_id}
        # 如果有上游父段，也排除它，避免共享端点误判。
        if vj.parent_id >= 0:
            exclude_vp.add(vj.parent_id)
        # vs 要排除原 vj 和旧孩子，因为它们在分叉点附近共享连接关系。
        exclude_vs = {vj_id, *vj.children}
        # vnew 只需要排除父段自身。
        exclude_vnew = {vj_id}

        # Exact finite-radius collision is checked after the temporary topology
        # has been installed and Murray radii are known.  A pre-insertion
        # centreline-intersection test is not authoritative in 3D and can
        # wrongly reject legal blended junctions.
        # 新终端分支还要满足最小间距约束。
        if self._segment_too_close_to_tree(vnew_xp, vnew_xd, exclude_vnew):
            return False

        # 所有局部约束都通过。
        return True

    def _passes_distal_constraints(self, vj_id: int, xnew: np.ndarray) -> bool:
        """
        check if a new branch can be attached to the distal end of the parent vessel `vj`.
        """
        vj = self.vessels[vj_id]
        if vj.branching_mode != "distal":
            return False
        
        # compute the length of the new segment. length must be greater than the minimum segment length.
        xb      = vj.x_d
        lnew    = _distance(xnew, xb)
        if lnew < 1e-6:
            return False
        if lnew < self.cfg.min_segment_length_um:
            return False
        
        # avoid creating a long spike: if the new terminal branch is too long, 
        # it usually indicates that the terminal point is too far from the tree, 
        # which is costly to connect and unnatural in morphology.
        if self.cfg.max_terminal_segment_um > 0 and lnew > self.cfg.max_terminal_segment_um:
            return False
        
        # check the bifurcation domain and perfusion domain constraints.
        if not (self._bifurcation_domain_ok(vj, xb) and self._in_perfusion_domain(xnew)):
            return False
        
        # new branch should not in banned region
        if not self._segment_domain_ok(xb, xnew, self.cfg.default_vessel_role):
            return False
        
        parent_dir  = vj.x_d - vj.x_p
        new_dir     = xnew - xb
        parent_bend = _angle_deg(parent_dir, new_dir)

        # check the bend angle constraints: the new branch should not bend too much relative to the parent vessel.
        if parent_bend > self.cfg.max_daughter_bend_deg:
            return False
        if parent_bend < self.cfg.min_new_branch_angle_deg:
            return False
        for c_id in vj.children:
            child_dir = self.vessels[c_id].x_d - self.vessels[c_id].x_p
            opening = _angle_deg(child_dir, new_dir)
            if opening < self.cfg.min_daughter_angle_deg or opening > self.cfg.max_opening_angle_deg:
                return False
            
        exclude = {vj_id}
        exclude.update(vj.children)
        # new branch should not be too close to the existing tree.
        if self._segment_too_close_to_tree(xb, xnew, exclude):
            return False
        
        return True

    def _aspect_ratios_ok(self, vessel_ids: Optional[Iterable[int]] = None) -> bool:
        """检查长径比约束：血管不能短到和半径同量级。"""
        # 默认检查全树，也可以只检查某些受影响 vessel。
        ids = vessel_ids if vessel_ids is not None else range(len(self.vessels))
        for vid in ids:
            # 忽略非法 id。
            if vid < 0 or vid >= len(self.vessels):
                continue
            v = self.vessels[vid]
            # 半径必须为正。
            if v.radius <= 0.0:
                return False
            # 长径比太小表示血管段短粗，几何上不合理。
            if v.length() / v.radius <= self.cfg.aspect_ratio_min:
                return False
        # 所有检查对象都通过。
        return True

    def _junction_symmetry_ok(self, vessel_ids: Optional[Iterable[int]] = None) -> bool:
        """
        avoid the situation where one child branch is much smaller than the other, which would indicate a severe asymmetry in the junction.
        """
        # delta <= 0 表示关闭该约束。
        if self.cfg.delta <= 0.0:
            return True
        # 默认检查全树。
        ids = vessel_ids if vessel_ids is not None else range(len(self.vessels))
        for vid in ids:
            # 忽略非法 id。
            if vid < 0 or vid >= len(self.vessels):
                continue
            # 只对至少两个孩子的 junction 检查半径对称性。
            children = self.vessels[vid].children
            if len(children) < 2:
                continue
            # 两两比较子分支半径比例。
            for i in range(len(children)):
                for j in range(i + 1, len(children)):
                    r1 = self.vessels[children[i]].radius
                    r2 = self.vessels[children[j]].radius
                    # 小半径/大半径太小，说明分支严重不对称。
                    if min(r1, r2) / max(r1, r2) <= self.cfg.delta:
                        return False
        # 所有 junction 对称性通过。
        return True

    def _murray_law_ok(self) -> bool:
        """Check Murray's law independently from optional daughter symmetry."""
        tolerance = 1.0e-10
        for vessel in self.vessels:
            if not vessel.children:
                continue
            if not np.isfinite(vessel.murray_residual):
                return False
            if vessel.murray_residual > tolerance:
                return False
        return True

    @staticmethod
    def _share_endpoint(first, second, tolerance=1.0e-8) -> bool:
        """Return true for legally connected segments that meet at one endpoint."""
        return any(
            np.linalg.norm(a - b) <= tolerance
            for a in (first.x_p, first.x_d)
            for b in (second.x_p, second.x_d)
        )

    def _topologically_local(self, first, second) -> bool:
        """
        Allow the natural blended junction region for segments at graph distance
        at most two.  Straight cylinders around a shared junction otherwise
        appear to overlap even though they form one continuous vascular lumen.
        """
        if self._share_endpoint(first, second):
            return True
        if first.parent_id >= 0 and first.parent_id == second.parent_id:
            return True
        first_parent = (
            self.vessels[first.parent_id] if first.parent_id >= 0 else None
        )
        second_parent = (
            self.vessels[second.parent_id] if second.parent_id >= 0 else None
        )
        return bool(
            (first_parent is not None and first_parent.parent_id == second.vid)
            or (second_parent is not None and second_parent.parent_id == first.vid)
        )

    def _radius_affected_ids(
        self, vessel_ids: Optional[Iterable[int]]
    ) -> set[int]:
        """Expand locally changed segments with their upstream Murray ancestors."""
        if vessel_ids is None:
            return set(range(len(self.vessels)))
        affected = {
            int(vid)
            for vid in vessel_ids
            if 0 <= int(vid) < len(self.vessels)
        }
        stack = list(affected)
        while stack:
            parent_id = self.vessels[stack.pop()].parent_id
            if parent_id >= 0 and parent_id not in affected:
                affected.add(parent_id)
                stack.append(parent_id)
        return affected

    def _tube_clearances_ok(
        self, vessel_ids: Optional[Iterable[int]] = None
    ) -> bool:
        """Check non-adjacent finite-radius tubes, not only centrelines."""
        clearance = float(self.cfg.min_vessel_clearance_um)
        affected = self._radius_affected_ids(vessel_ids)
        checked = set()
        for first_index in affected:
            first = self.vessels[first_index]
            for second_index, second in enumerate(self.vessels):
                if first_index == second_index:
                    continue
                pair = (
                    min(first_index, second_index),
                    max(first_index, second_index),
                )
                if pair in checked:
                    continue
                checked.add(pair)
                if self._topologically_local(first, second):
                    continue
                distance = float(
                    segment_distances_3d(
                        first.x_p,
                        first.x_d,
                        np.asarray([second.x_p], dtype=float),
                        np.asarray([second.x_d], dtype=float),
                    )[0]
                )
                required = first.radius + second.radius + clearance
                if distance + 1.0e-9 < required:
                    return False
        return True

    def _tube_boundaries_ok(
        self, vessel_ids: Optional[Iterable[int]] = None
    ) -> bool:
        """Keep each complete circular tube inside the rectangular domain."""
        surface_clearance = float(self.cfg.boundary_surface_clearance_um)
        for vid in self._radius_affected_ids(vessel_ids):
            vessel = self.vessels[vid]
            margin = vessel.radius + surface_clearance
            # The authoritative root begins on the inlet face by construction.
            if vessel.parent_id >= 0 and not self._in_bounds(vessel.x_p, margin):
                return False
            if not self._in_bounds(vessel.x_d, margin):
                return False
            excluded_regions = list(self.cfg.avascular_regions)
            if vessel.role not in {"transport", "perforator"}:
                excluded_regions.extend(self.cfg.carriage_regions)
            if excluded_regions:
                sample_count = max(
                    9,
                    int(np.ceil(vessel.length() / max(vessel.radius, 1.0))) + 1,
                )
                for t in np.linspace(0.0, 1.0, sample_count):
                    point = (1.0 - t) * vessel.x_p + t * vessel.x_d
                    if any(
                        region.signed_distance(point) + 1.0e-9 < vessel.radius
                        for region in excluded_regions
                    ):
                        return False
        return True

    def _post_radius_constraints_ok(
        self,
        angle_vessel_ids: Optional[Iterable[int]] = None,
    ) -> bool:
        """
        Check all radius-dependent constraints after the direct Murray update.
        """
        # 长径比检查默认全树，因为半径更新可能影响许多祖先段。
        return (
            self._aspect_ratios_ok()
            and self._murray_law_ok()
            and self._junction_symmetry_ok()
            and self._junction_geometry_ok(angle_vessel_ids)
            and self._tube_boundaries_ok(angle_vessel_ids)
            and self._tube_clearances_ok(angle_vessel_ids)
        )
