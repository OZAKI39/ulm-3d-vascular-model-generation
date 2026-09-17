"""
血管树生长主流程。

这个模块回答“树怎么一步一步长出来”：
1. 初始化第一条 root 血管；
2. 随机采样一个尚未被血管覆盖的新终端点；
3. 在现有树中寻找适合连接它的父血管；
4. 临时插入分叉、更新半径、检查约束；
5. 如果合法就保留，否则回滚并尝试下一个候选。

为什么要“临时插入再回滚”：
DCCO 的很多约束依赖插入后的全局状态，例如半径更新后是否碰撞、
分叉角是否仍合法、体积代价是否变小。只看插入前的几何距离不够可靠，
所以必须先模拟一次真实插入，再决定是否接受。
"""

import math
import numpy as np
from ..geometry.geometry import _distance
from ..core.models import Vessel


class GrowthMixin:
    """
    Responsible for initializing the DCCO tree, sampling new terminal points, and growing the tree step by step.
    """
    def _terminal_outflow(self):
        """
        Return a placeholder; _apply_radii allocates equal terminal flow globally.
        """
        return 0.0

    def _initialize_root(self):
        """
        initialize the first root vessel of the entire tree.
        """

        # get the proximal coordinates of the root vessel from the configuration. 
        # If they are NaN, sample a random point on the boundary.
        root_x = self.cfg.root_x_um
        root_y = self.cfg.root_y_um
        root_z = self.cfg.root_z_um
        if self.cfg.geometry_mode == "planar_2d":
            root_y = self.cfg.plane_y_um
            if math.isnan(root_x) or math.isnan(root_z):
                side = self.cfg.cube_um
                margin = 0.2 * side
                edge = int(self.rng.integers(0, 4))
                coordinate = float(self.rng.uniform(margin, side - margin))
                if edge == 0:
                    root_x, root_z = 0.0, coordinate
                elif edge == 1:
                    root_x, root_z = side, coordinate
                elif edge == 2:
                    root_x, root_z = coordinate, 0.0
                else:
                    root_x, root_z = coordinate, side
        elif math.isnan(root_x) or math.isnan(root_y) or math.isnan(root_z):
            lengths = self._domain_lengths()
            margins = 0.2 * lengths
            point = np.asarray(
                [
                    self.rng.uniform(margins[0], lengths[0] - margins[0]),
                    self.rng.uniform(margins[1], lengths[1] - margins[1]),
                    self.rng.uniform(margins[2], lengths[2] - margins[2]),
                ],
                dtype=float,
            )
            face = int(self.rng.integers(0, 6))
            axis = face // 2
            point[axis] = 0.0 if face % 2 == 0 else lengths[axis]
            root_x, root_y, root_z = point
        x_p = np.array([root_x, root_y, root_z], dtype=float)

        # x_d_best is the best distal point found for the root vessel. 
        # It is initialized to None and will be updated if a legal distal point is found.
        x_d_best = None

        # randomly choose a distal point for the root vessel within the domain, 
        # trying up to 4000 times to find a legal point.
        for _ in range(4000):
            candidate = self._random_point()

            # distal point must be within the perfusion domain.
            if not self._in_perfusion_domain(candidate):
                continue

            # length of the root segment must be within the configured minimum and maximum lengths.
            # TODO: more reliable way to decide the length
            length = _distance(candidate, x_p)
            if not (self.cfg.min_root_segment_um <= length <= self.cfg.max_root_segment_um):
                continue

            # root segment must be within the allowed domain for the default vessel role.
            if not self._segment_domain_ok(x_p, candidate, self.cfg.default_vessel_role):
                continue
            
            # Segment aspect ratio uses half-width in 2-D and radius in 3-D.
            # TODO: the parameter is from formula (5) in paper, show the derivation
            if (
                length / self.haemodynamics.root_size_um
                > self.cfg.aspect_ratio_min
            ):
                x_d_best = candidate
                break
        
        if x_d_best is None:
            raise RuntimeError(
                "Could not initialize the root vessel because no valid distal "
                "point was found after 4000 attempts."
            )

        # create the root vessel with the chosen proximal and distal points, and add it to the list of vessels.
        root = Vessel(
            vid=0,
            parent_id=-1,
            x_p=x_p,
            x_d=x_d_best,
            radius=self.haemodynamics.root_size_um,
            branching_mode=self.cfg.default_branching_mode,
            role=self.cfg.default_vessel_role,
            prescribed_outflow=0.0,
        )

        # add the root vessel to the list of vessels, invalidate the tree cache, and apply the radii to the vessels.
        self.vessels.append(root)
        self._invalidate_tree_cache()
        self._apply_radii()

    def _add_terminal(self):
        """
        try to add a new terminal point to the existing tree.
        """
        for _ in range(40):
            # sample a new terminal point that is not too close to existing vessels and is within the perfusion domain.
            x_new = self._sample_terminal_point()
            if x_new is None:
                continue

            # insert the new terminal point into the existing tree by finding the best parent vessel and bifurcation point.
            success = self._try_insert_terminal(x_new)
            if success:
                self._fr_runtime = 1.0
                return True
            
        return False

    def _sample_terminal_point(self):
        """
        Sample a new terminal point that is not too close to existing vessels, and is within the perfusion domain.
        It should not be too far away from the existing tree, otherwise connecting it would create a long spiky vessel.
        TODO: find the minimum distance threshold from the existing tree
        """
        n_fail  = 0
        upper   = self.cfg.max_sample_distance_um
        
        # try maximum 600 times to find a valid terminal point 
        for _ in range(600):
            p = self._random_point()

            # check if the sampled point is within the perfusion domain. 
            next_terminal_count = max(self._terminal_count() + 1, 1)
            predicted_terminal_radius = self.haemodynamics.root_size_um * (
                1.0 / next_terminal_count
            ) ** (1.0 / self.haemodynamics.murray_gamma)
            boundary_margin = (
                predicted_terminal_radius
                + float(self.cfg.boundary_surface_clearance_um)
            )
            valid_domain = (
                self._in_perfusion_domain(p)
                and self._in_bounds(p, boundary_margin)
            )

            # check the closest distance from the sampled point to the existing vessels.
            d = self._closest_distance(p) if valid_domain else math.inf

            # new point distance to the existing vessels must be greater than the minimum distance threshold,
            # and it must not exceed the configured maximum connection distance.
            if valid_domain and d > self._lmin() and (upper <= 0 or d <= upper):
                return p
            
            n_fail += 1
            if n_fail >= self.cfg.n_fail:
                # gradually weakens the minimum-distance constraint to make terminal generation easier in crowded regions
                self._fr_runtime *= self.cfg.fr
                n_fail = 0
                if self._fr_runtime < 1e-6:
                    return None
        
        return None

    def _try_insert_terminal(self, x_new):
        """
        choose the best parent vessel and bifurcation point to connect the new terminal point x_new to the existing tree.
        """
        # Record 2-D lumen area or 3-D lumen volume before insertion.
        baseline_lumen_measure = self._lumen_measure_cost()

        # decide the new terminal's tissue demand once
        terminal_outflow = self._terminal_outflow()

        # find the candidate parent vessels in the neighborhood of the new terminal point.
        neighbourhood = self._neighborhood(x_new)

        best_cost   = math.inf
        best_choice = None

        # try inserting each candidate bifurcation point for each neighboring parent vessel and calculate the cost.
        for vj_id in neighbourhood:
            vj = self.vessels[vj_id]
            for xb in self._candidate_bifurcation_points(vj, x_new):
                cost = self._trial_bifurcation_cost(
                    vj_id,
                    xb,
                    x_new,
                    baseline_lumen_measure,
                    terminal_outflow,
                )
                if cost is None:
                    continue

                # save the best candidate with the lowest cost that satisfies all constraints.
                if cost < best_cost:
                    best_cost = cost
                    best_choice = (vj_id, xb)

        # no valid candidate found, insertion fails.
        if best_choice is None:
            return False

        vj_id, xb = best_choice
        if self.vessels[vj_id].branching_mode == "distal":
            # check if the new point satisfies the distal constraints before attaching it.
            if not self._passes_distal_constraints(vj_id, x_new):
                return False
            
            # update the tree structure by attaching the new terminal point to the distal end of the parent vessel.
            vnew_id, undo = self._attach_distal(vj_id, x_new, terminal_outflow)
            self._apply_radii()

            # check if the updated tree satisfies the post-radius constraints after the distal attachment.
            if not self._post_radius_constraints_ok({vj_id}):
                self._undo_distal_attach(vj_id, vnew_id, undo)
                self._apply_radii()
                return False
        else:
            if not self._passes_local_constraints(vj_id, xb, x_new):
                return False
            
            vs_id, vnew_id, undo = self._split_vessel(vj_id, xb, x_new, terminal_outflow)

            affected = {vj_id, vs_id}
            if undo["old_children"]:
                affected.update(undo["old_children"])
            if self.vessels[vj_id].parent_id >= 0:
                affected.add(self.vessels[vj_id].parent_id)

            self._apply_radii()

            if not self._post_radius_constraints_ok(affected):
                self._undo_split(vj_id, vs_id, vnew_id, undo)
                self._apply_radii()
                return False

        self._invalidate_tree_cache()
        return True

    def _trial_bifurcation_cost(
        self,
        vj_id,
        xb,
        x_new,
        baseline_lumen_measure,
        terminal_outflow,
    ):
        """
        compute the cost of a trial bifurcation.
        return None if the candidate violates any constraints; otherwise, return the objective function increment.
        """
        if self.cfg.cost_function != "sprouting":
            raise ValueError("cost_function must be 'sprouting'.")

        # vj is the current candidate parent vessel.
        vj = self.vessels[vj_id]

        if vj.branching_mode == "distal":
            # Prescribed target outflows occur only at leaves.
            if vj.prescribed_outflow > 0.0:
                return None

            # do the local constraints check first, which does not require modifying the tree structure.
            if not self._passes_distal_constraints(vj_id, x_new):
                return None

            branch_point = vj.x_d.copy()

            # modify the tree structure temporarily 
            vnew_id, undo = self._attach_distal(vj_id, x_new, terminal_outflow)
            try:
                try:
                    self._apply_radii()
                except RuntimeError:
                    return None
                if not self._post_radius_constraints_ok({vj_id}):
                    return None
                trial_lumen_measure = self._lumen_measure_cost()
            finally:
                self._undo_distal_attach(vj_id, vnew_id, undo)
                self._apply_radii()
        else:
            # check local constraints for the split mode
            if not self._passes_local_constraints(vj_id, xb, x_new):
                return None

            branch_point = xb
            vs_id, vnew_id, undo = self._split_vessel(vj_id, xb, x_new, terminal_outflow)
            try:
                affected = {vj_id, vs_id}
                if undo["old_children"]:
                    affected.update(undo["old_children"])
                if vj.parent_id >= 0:
                    affected.add(vj.parent_id)
                try:
                    self._apply_radii()
                except RuntimeError:
                    return None
                if not self._post_radius_constraints_ok(affected):
                    return None
                trial_lumen_measure = self._lumen_measure_cost()
            finally:
                self._undo_split(vj_id, vs_id, vnew_id, undo)
                self._apply_radii()

        lumen_measure_delta = trial_lumen_measure - baseline_lumen_measure
        
        vref = max(
            self.haemodynamics.reference_lumen_measure(self._lc()),
            1e-12,
        )
        lref = max(self._lc(), 1e-12)                   
        rref = max(self.haemodynamics.root_size_um, 1e-12)
        lnew = _distance(x_new, branch_point)

        return (
            self.cfg.sprout_cv * lumen_measure_delta / vref
            + self.cfg.sprout_cp * vj.radius / rref         # penalty for sprouting from a thick parent vessel
            + self.cfg.sprout_cd * (lnew / lref) ** 2       # penalty for new branches that are too long
        )

    def generate(self):
        """
        This is for generating a complete vascular tree starting from an empty tree.
        Not for staged growth, which is handled by `generate_staged_tree()`.
        """
        self._initialize_root()
        consecutive_failures = 0

        while self._terminal_count() < self.cfg.n_terminals:
            if self._add_terminal():
                consecutive_failures = 0
                continue

            consecutive_failures    += 1
            self._fr_runtime        = 1.0

            if consecutive_failures >= self.cfg.max_terminal_insert_failures:
                break
        if self._terminal_count() != int(self.cfg.n_terminals):
            raise RuntimeError(
                "Vascular growth stopped before the configured terminal count "
                f"was reached: generated={self._terminal_count()}, "
                f"target={int(self.cfg.n_terminals)}."
            )
        self._apply_radii()
        return self.vessels

    def grow_additional_terminals(self, n_terminals):
        """
        Add designated number of terminals to the existing tree. 
        This is used for staged growth: 
            each stage can continue growing on the previous tree, or grow independently on a Tinit subtree and then merge.
        """
        added = 0
        if not self.vessels and n_terminals > 0:
            self._initialize_root()
            added += 1

        consecutive_failures = 0
        while added < n_terminals:
            if self._add_terminal():
                added += 1
                consecutive_failures = 0
            else:
                consecutive_failures += 1
                self._fr_runtime = 1.0

            if consecutive_failures >= self.cfg.max_terminal_insert_failures:
                break

        if added != int(n_terminals):
            raise RuntimeError(
                "Staged vascular growth stopped before its requested number of "
                f"terminals was added: added={added}, requested={int(n_terminals)}."
            )
        self._apply_radii()
        return added
