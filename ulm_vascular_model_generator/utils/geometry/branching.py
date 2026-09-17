"""
分叉操作模块。

这个模块只负责“怎样把新血管接到旧血管上”。在 DCCO 中，新增终端时
通常不是简单把新段挂到父血管末端，而是把父血管从某个分叉点 `xb`
切开，形成：
- `vp`：父血管的近端部分，继续沿用原来的 vessel id；
- `vs`：父血管的远端延续段，继承旧远端和旧子树；
- `vnew`：从分叉点连接到新终端点的新分支。

为什么要保留 undo 信息：
候选分叉点要经过试插入和约束检查。很多候选会失败，因此每次修改树
之前都保存旧状态，失败后能精确回滚。
"""

from __future__ import annotations
import numpy as np
from .geometry import _angles_deg, _distance
from ..core.models import Vessel


class BranchingMixin:
    """集中处理血管分裂、远端挂接和候选分叉点生成。"""

    def _new_branching_mode(self) -> str:
        """新生成分支使用 YAML 中的默认分叉模式。"""
        # 新增 vnew 的 branching_mode 由全局默认配置决定。
        return self.cfg.default_branching_mode

    def _split_vessel(self, vj_id: int, xb: np.ndarray, xnew: np.ndarray,
                       prescribed_outflow: float = 0.0) -> tuple[int, int, dict]:
        """
        把一条父血管 `vj` 在 `xb` 处分裂，并接入新终端 `xnew`。

        返回 `(vs_id, vnew_id, undo_state)`。调用方可以用 undo_state 回滚。
        这里叫 “Permanently split”，但在试算阶段也会调用；是否最终保留
        由调用方决定。
        """
        # vj 是即将被切开的父血管。
        vj = self.vessels[vj_id]
        # 保存旧状态。注意不仅要保存几何和子节点，还要保存血流相关字段，
        # 因为试插入会调用 `_apply_radii()` 改写这些字段。
        undo = {
            "old_x_d": vj.x_d.copy(),
            "old_children": list(vj.children),
            "old_radius": vj.radius,
            "old_flow_rate": vj.flow_rate,
            "old_mean_velocity": vj.mean_velocity,
            "old_flow_conservation_residual": vj.flow_conservation_residual,
            "old_murray_residual": vj.murray_residual,
            "old_prescribed_outflow": vj.prescribed_outflow,
        }

        # vs 会追加到列表末尾，因此新 id 是当前长度。
        vs_id = len(self.vessels)
        # vs 是旧父血管的远端延续段。它继承旧父血管的远端、旧孩子、
        # branching_mode、role，以及可选主干标记。
        vs = Vessel(
            vid=vs_id,
            parent_id=vj_id,
            x_p=xb.copy(),
            x_d=undo["old_x_d"].copy(),
            radius=vj.radius,
            children=list(undo["old_children"]),
            branching_mode=vj.branching_mode,
            role=vj.role,
            is_main_trunk=vj.is_main_trunk,
            prescribed_outflow=vj.prescribed_outflow,
        )
        # 把远端延续段加入树。
        self.vessels.append(vs)
        # 旧孩子原来挂在 vj 上；分裂后它们应该挂到远端延续段 vs 上。
        for c in undo["old_children"]:
            self.vessels[c].parent_id = vs_id

        # vnew 也追加到列表末尾。
        vnew_id = len(self.vessels)
        # vnew 是真正新增的终端分支。
        vnew = Vessel(
            vid=vnew_id,
            parent_id=vj_id,
            x_p=xb.copy(),
            x_d=xnew.copy(),
            radius=self.haemodynamics.minimum_size_um,
            branching_mode=self._new_branching_mode(),
            role=self.cfg.default_vessel_role,
            is_main_trunk=False,
            prescribed_outflow=max(float(prescribed_outflow), 0.0),
        )
        # 把新终端分支加入树。
        self.vessels.append(vnew)

        # vj 本身变成近端段 vp：远端点变为 xb，孩子变为 [vs, vnew]。
        vj.x_d = xb.copy()
        vj.children = [vs_id, vnew_id]
        # 如果 vj 原来携带 prescribed_outflow（例如主干 outlet），这个流量
        # 应该跟随旧远端移动到 vs，而不是留在近端段上重复计入。
        vj.prescribed_outflow = 0.0

        # 返回两个新段 id 和回滚所需的旧状态。
        return vs_id, vnew_id, undo

    def _undo_split(self, vj_id: int, vs_id: int, vnew_id: int, undo: dict) -> None:
        """撤销最近一次 `_split_vessel`。"""
        # 找回被切开的原父血管。
        vj = self.vessels[vj_id]
        # 恢复父血管原始几何、子节点和血流字段。
        vj.x_d = undo["old_x_d"]
        vj.children = undo["old_children"]
        vj.radius = undo["old_radius"]
        vj.flow_rate = undo["old_flow_rate"]
        vj.mean_velocity = undo["old_mean_velocity"]
        vj.flow_conservation_residual = undo["old_flow_conservation_residual"]
        vj.murray_residual = undo["old_murray_residual"]
        vj.prescribed_outflow = undo["old_prescribed_outflow"]
        # 旧孩子的父节点重新指回 vj。
        for c in undo["old_children"]:
            self.vessels[c].parent_id = vj_id
        # `_split_vessel` 保证 vs 和 vnew 是最后 append 的两个对象，
        # 所以回滚时可以从列表末尾 pop，避免移动其它 vessel 的 id。
        # 先删除最后追加的 vnew。
        if vnew_id == len(self.vessels) - 1:
            self.vessels.pop()
        # 再删除倒数第二个追加的 vs。
        if vs_id == len(self.vessels) - 1:
            self.vessels.pop()

    def _attach_distal(self, vj_id: int, xnew: np.ndarray,
                       prescribed_outflow: float = 0.0) -> tuple[int, dict]:
        """
        distal 模式：不切开父血管，而是在父血管远端直接挂一条新分支。
        """
        vj = self.vessels[vj_id]
        undo = {"old_children": list(vj.children)}
        vnew_id = len(self.vessels)

        self.vessels.append(Vessel(
            vid=vnew_id,
            parent_id=vj_id,
            x_p=vj.x_d.copy(),
            x_d=xnew.copy(),
            radius=self.haemodynamics.minimum_size_um,
            branching_mode=self._new_branching_mode(),
            role=self.cfg.default_vessel_role,
            is_main_trunk=False,
            prescribed_outflow=max(float(prescribed_outflow), 0.0),
        ))

        vj.children.append(vnew_id)

        return vnew_id, undo

    def _undo_distal_attach(self, vj_id: int, vnew_id: int, undo: dict) -> None:
        """撤销 distal 挂接。"""
        # 找回父血管。
        vj = self.vessels[vj_id]
        # 恢复旧 children。
        vj.children = undo["old_children"]
        # 如果新分支仍在列表末尾，就直接 pop 删除。
        if vnew_id == len(self.vessels) - 1:
            self.vessels.pop()

    def _triangle_points(self, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray
                          ) -> list[np.ndarray]:
        """
        generate bifurcation points inside the triangle formed by p1, p2, p3.
        """
        # control the density of candidate points inside the triangle.
        n = max(self.cfg.delta_v, 3)

        i, j = np.meshgrid(np.arange(1, n), np.arange(1, n), indexing="ij")
        mask = (i + j) < n
        if not np.any(mask):
            return []
        
        ii = i[mask].astype(float)
        jj = j[mask].astype(float)
        kk = n - ii - jj

        bary = np.column_stack((ii, jj, kk)) / float(n)
        pts = bary[:, 0, None] * p1 + bary[:, 1, None] * p2 + bary[:, 2, None] * p3
        return [pts[row].copy() for row in range(len(pts))]

    def _line_points(self, p1: np.ndarray, p2: np.ndarray) -> list[np.ndarray]:
        """
        generate bifurcation points along the line segment from p1 to p2.
        suitable for fixed branching mode, where the bifurcation point is constrained to lie on the parent vessel segment.
        """
        n = max(self.cfg.delta_v, 3)

        t = (np.arange(1, n, dtype=float) / float(n))[:, None]
        pts = (1.0 - t) * p1 + t * p2
        return [pts[row].copy() for row in range(len(pts))]

    def _prefilter_bifurcation_points(
        self,
        vj: Vessel,
        xnew: np.ndarray,
        points: list[np.ndarray],
    ) -> list[np.ndarray]:
        """
        Filter the candidate bifurcation points based on length and angle constraints.
        Reduce the number of expensive trial insertions by eliminating obviously invalid points.
        """
        if not points:
            return []

        pts = np.asarray(points, dtype=float)

        vp_dir      = pts - vj.x_p
        vs_dir      = vj.x_d - pts
        vnew_dir    = xnew - pts

        lp      = np.linalg.norm(vp_dir, axis=1)
        ls      = np.linalg.norm(vs_dir, axis=1)
        lnew    = np.linalg.norm(vnew_dir, axis=1)

        vj_length = _distance(vj.x_d, vj.x_p)
        min_parent_split_segment = max(
            self.cfg.min_segment_length_um,
            self.cfg.min_split_fraction * vj_length,
        )

        # ========== the first round of filtering: length constraints. ==========
        # the new bifurcation point must be far enough from the parent vessel endpoints, 
        # and the new branch must be long enough.
        mask = (
            (lp >= min_parent_split_segment)
            & (ls >= min_parent_split_segment)
            & (lnew >= self.cfg.min_segment_length_um)
        )

        # limit the maximum length of the new branch if specified in the configuration.
        if self.cfg.max_terminal_segment_um > 0:
            mask &= lnew <= self.cfg.max_terminal_segment_um

        # no candidates pass the length filter, return early.
        if not np.any(mask):
            return []

        # ========== the second round of filtering: angle constraints. ==========
        bend_vs     = _angles_deg(vp_dir, vs_dir)
        bend_vnew   = _angles_deg(vp_dir, vnew_dir)
        opening     = _angles_deg(vs_dir, vnew_dir)
        
        # the new bifurcation point must satisfy the angle constraints for the parent vessel, the new branch, and the opening angle between the two branches.
        mask &= (
            (bend_vs <= self.cfg.max_daughter_bend_deg)
            & (bend_vnew <= self.cfg.max_daughter_bend_deg)
            & (bend_vnew >= self.cfg.min_new_branch_angle_deg)
            & (opening >= self.cfg.min_daughter_angle_deg)
            & (opening <= self.cfg.max_opening_angle_deg)
        )
        # no candidates pass the angle filter, return early.
        if not np.any(mask):
            return []

        # the directions of the existing children must also be consistent with the new distal continuation segment.
        for c_id in vj.children:
            child = self.vessels[c_id]
            child_dir = child.x_d - child.x_p
            mask &= _angles_deg(vs_dir, child_dir) <= self.cfg.max_daughter_bend_deg
            if not np.any(mask):
                return []

        if vj.parent_id >= 0:
            parent      = self.vessels[vj.parent_id]
            parent_dir  = parent.x_d - parent.x_p
            mask &= _angles_deg(vp_dir, parent_dir) <= self.cfg.max_daughter_bend_deg
            if not np.any(mask):
                return []

            for sibling_id in parent.children:
                if sibling_id == vj.vid:
                    continue
                sibling         = self.vessels[sibling_id]
                sibling_dir     = sibling.x_d - sibling.x_p
                sibling_opening = _angles_deg(vp_dir, sibling_dir)
                mask &= (
                    (sibling_opening >= self.cfg.min_daughter_angle_deg)
                    & (sibling_opening <= self.cfg.max_opening_angle_deg)
                )
                if not np.any(mask):
                    return []

        # apply the mask
        filtered = pts[mask]
        return [filtered[row].copy() for row in range(len(filtered))]

    def _candidate_bifurcation_points(self, vj: Vessel, xnew: np.ndarray) -> list[np.ndarray]:
        """
        return a list of candidate bifurcation points for a new branch at the parent vessel vj.
        """
        # non_branching does not allow any bifurcation points.
        if vj.branching_mode == "non_branching":
            return []
        # distal mode only allows bifurcation at the distal end of the parent vessel.
        if vj.branching_mode == "distal":
            return [vj.x_d.copy()]
        
        if vj.branching_mode == "fixed":
            # fixed mode only allows bifurcation along the parent vessel segment.
            points = self._line_points(vj.x_p, vj.x_d)
        else:
            # versatile mode allows bifurcation points to move inside the triangle formed by xnew, vj.x_p, and vj.x_d.
            points = self._triangle_points(xnew, vj.x_p, vj.x_d)
        
        return self._prefilter_bifurcation_points(vj, xnew, points)
