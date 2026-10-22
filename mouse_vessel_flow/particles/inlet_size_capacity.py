"""Certified inlet clearance bounds and exact flux proposal acceleration.

Only the real open inlet triangles and solid WALL enter this calculation.
Distance to a closed set is 1-Lipschitz. A centroid distance +/- the furthest
vertex distance therefore bounds every point of a triangle. Refinement can
discard only triangles proved impossible; it never accepts a birth by itself.
"""
import heapq
import math
import numpy as np
from .inlet_flux import InletFluxSampler
from .particle82a_geometry import maximum_handoff_radius, lower_gap


class InletClearanceTree:
    def __init__(self, sampler, wall):
        self.sampler, self.wall = sampler, wall
        self.nodes = []
        self.ro = float(wall.roundoff_m)
        self.distance_evaluations = 0
        self.roots = [self.add(x, q, i, 0, (k,)) for k, (i, x, q, _) in enumerate(sampler.pieces)]

    def add(self, xyz, q, original_triangle, depth, path):
        center = xyz.mean(axis=0)
        distance = self.wall.nearest_center_triangle(center)[1]
        self.distance_evaluations += 1
        rho = float(np.linalg.norm(xyz-center, axis=1).max())
        # Existing exact positive P1 flux definition, including clipped pieces.
        weight = float(np.linalg.norm(np.cross(xyz[1]-xyz[0], xyz[2]-xyz[0]))/2*q.mean())
        node = dict(xyz=xyz, q=q, original_triangle=int(original_triangle),
                    center=center, distance=float(distance),
                    lower=max(0., distance-rho-self.ro), upper=distance+rho+self.ro,
                    weight=weight, children=None, depth=depth, path=path)
        self.nodes.append(node)
        return len(self.nodes)-1

    def split(self, key):
        node = self.nodes[key]
        if node['children'] is None:
            x, q = node['xyz'], node['q']
            pairs = [(0, 1), (1, 2), (2, 0)]
            i, j = max(pairs, key=lambda ij: np.linalg.norm(x[ij[0]]-x[ij[1]]))
            k = 3-i-j
            m, qm = (x[i]+x[j])/2, (q[i]+q[j])/2
            node['children'] = [self.add(np.array([x[i], m, x[k]]), np.array([q[i], qm, q[k]]), node['original_triangle'], node['depth']+1, node['path']+(0,)),
                                self.add(np.array([m, x[j], x[k]]), np.array([qm, q[j], q[k]]), node['original_triangle'], node['depth']+1, node['path']+(1,))]
        return node['children']

    def capacity(self, tolerance_m=1e-12, max_refinements=200000):
        """Global maximum bracket over the whole positive-flux cap, not a raster maximum."""
        heap = [(-self.nodes[i]['upper'], i) for i in self.roots]
        heapq.heapify(heap)
        witness = max(self.roots, key=lambda i: self.nodes[i]['distance'])
        best = self.nodes[witness]['distance']
        refinements = 0
        while -heap[0][0]-best > tolerance_m:
            _, key = heapq.heappop(heap)
            for child in self.split(key):
                n = self.nodes[child]
                if n['distance'] > best:
                    best, witness = n['distance'], child
                heapq.heappush(heap, (-n['upper'], child))
            refinements += 1
            if refinements >= max_refinements:
                raise RuntimeError('GLOBAL_CAPACITY_NUMERICAL_PROGRESS_FAILURE')
        n = self.nodes[witness]
        # Include only the existing admission roundoff allowance; no new gap.
        lo = float(maximum_handoff_radius(best+self.ro))
        hi = float(maximum_handoff_radius(-heap[0][0]+self.ro))
        return dict(schema='CERTIFIED_INLET_CAPACITY_V1',
            D_source_max_m=4e-6, D_geometry_max_m=2*lo,
            D_geometry_max_bracket_m=[2*lo, 2*hi], radius_max_bracket_m=[lo, hi],
            maximum_globally_passable_inlet_radius_m=lo,
            maximum_globally_passable_inlet_diameter_m=2*lo,
            pure_wall_clearance_diameter_bracket_m=[2*best, -2*heap[0][0]],
            witness_center_m=n['center'].tolist(), witness_triangle=n['original_triangle'],
            witness_inward_velocity_m_s=float(n['q'].mean()),
            wall_clearance_bracket_m=[best, -heap[0][0]],
            tolerance_m=tolerance_m, roundoff_m=self.ro, refinements=refinements,
            distance_evaluations=self.distance_evaluations,
            method='GLOBAL_TRIANGLE_BRANCH_AND_BOUND_OF_EXACT_WALL_DISTANCE; 1_LIPSCHITZ',
            scope='OPEN_INLET_BIRTH_ADMISSION_WITH_EXISTING_HANDOFF; NOT_DOWNSTREAM_VESSEL_PASSABILITY',
            positive_flux_cap_covers_entire_cap=bool(np.all(self.sampler.q >= 0) and np.all(self.sampler.q.sum(axis=1)>0)),
            upper_bound_is_not_treated_as_a_feasible_witness=True)

    def feasible_proposal(self, radius, target_certified_fraction=.125, max_refinements=200000):
        """Return exact original flux conditioned on a certified superset of feasibility.

        The final authoritative admission rejection then conditions on feasibility
        itself. The retained set contains ALL feasible points of positive flux.
        """
        required = float(radius+lower_gap(radius)-self.ro)
        full, pending = [], []
        def insert(key):
            n = self.nodes[key]
            if n['upper'] < required:
                return
            if n['lower'] >= required:
                full.append(key)
            else:
                heapq.heappush(pending, (-n['weight'], n['path'], key))
        for key in self.roots:
            insert(key)
        refinements = 0
        while True:
            good = math.fsum(self.nodes[k]['weight'] for k in full)
            ambiguous = math.fsum(-w for w, _, _ in pending)
            if good > 0 and good/(good+ambiguous) >= target_certified_fraction:
                break
            if not pending or refinements >= max_refinements:
                raise RuntimeError('FEASIBLE_FLUX_BOUND_NUMERICAL_PROGRESS_FAILURE')
            _, _, key = heapq.heappop(pending)
            for child in self.split(key):
                insert(child)
            refinements += 1
        keys = sorted(full+[k for _, _, k in pending], key=lambda k: self.nodes[k]['path'])
        nodes = [self.nodes[k] for k in keys]
        proposal = InletFluxSampler(np.array([n['xyz'] for n in nodes]), np.array([n['q'] for n in nodes]))
        fraction = good/proposal.Q_m3_s
        # Extra factor 4 is execution slack, not a geometric safety margin.
        # P(no success) <= 2^-64 at the certified geometric acceptance bound.
        guard = int(math.ceil(4*math.log(2**-64)/math.log1p(-min(fraction, 1-1e-15))))
        record = dict(feasible_flux_fraction_lower=good/self.sampler.Q_m3_s,
                      feasible_flux_fraction_upper=proposal.Q_m3_s/self.sampler.Q_m3_s,
                      accelerated_acceptance_lower=fraction, retained_triangle_count=len(keys),
                      position_guard=guard, guard_tail_probability_target=2**-64,
                      guard_execution_slack=4, refinement_count=refinements,
                      acceleration='EXACT_P1_FLUX_ON_CERTIFIED_SUPERSET; FINAL_ORIGINAL_ADMISSION_REJECTION')
        return proposal, np.array([n['original_triangle'] for n in nodes]), record


def compute_inlet_global_size_capacity(sampler, wall, **kwargs):
    return InletClearanceTree(sampler, wall).capacity(**kwargs)
