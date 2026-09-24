"""Memoize repeated identical geometry queries without changing P6.5 bytecode.

Uses P6's existing private dependency-binding mechanism. No module global or
upstream class is patched. Cache keys contain exact float64 positions/radii,
and cached geometry results are only read by the unchanged solver functions.
"""
from functools import lru_cache
from copy import copy
from types import MethodType
import numpy as np
from .particle_shapes import Sphere,roundoff_length
from .wall_gap import wall_gap
from .nearfield_handoff import exact_interactions,wall_handoff_certificate
from .physical_time_refinement import swept_clearance_certificate
from .particle4_motion import initial_world
from .particle65_motion import Particle65Stepper,v1_trial,assemble_v1
from .particle6_stepper import bind_query_dependency


def cached_step_to(stepper,*,global_motion_bound=False):
    original_wall=stepper.wall
    wall=original_wall
    if wall is not None:
        # Per-stepper memoization of EXACT BVH requests. The underlying triangles,
        # locator, normals and all geometric arithmetic remain unchanged.
        wall=copy(original_wall)
        @lru_cache(maxsize=128)
        def candidates_at(center,radius):
            ids=original_wall.candidates(center,radius)
            ids.flags.writeable=False
            return ids
        @lru_cache(maxsize=128)
        def nearest_at(center):return original_wall.nearest_center_triangle(center)
        wall.candidates=lambda center,radius:candidates_at(tuple(center),float(radius))
        wall.nearest_center_triangle=lambda center:nearest_at(tuple(center))
        stepper.wall=wall
    @lru_cache(maxsize=64)
    def gap_at(center,radius,inside):
        return wall_gap(Sphere(center,radius),wall,inside_lumen=inside)
    def cached_gap(shape,query_wall,*,inside_lumen=None):
        if query_wall is not wall or type(shape) is not Sphere:
            return wall_gap(shape,query_wall,inside_lumen=inside_lumen)
        return gap_at(tuple(shape.center_m),shape.radius_m,inside_lumen)
    def global_clearance(start,end,query_wall,lower):
        if query_wall is not wall or type(start) is not Sphere or type(end) is not Sphere or start.radius_m!=end.radius_m:return None
        # Distance to the union of original finite triangles is 1-Lipschitz.
        # This is the same surface-motion bound used by P3/P6.5, applied to
        # the global nearest gap before enumerating individual swept faces.
        gap=cached_gap(start,wall).gap_m
        bound=gap-float(np.linalg.norm(end.center_m-start.center_m))-lower
        pad=max(wall.roundoff_m,roundoff_length(start.center_m,end.center_m,start.radius_m))
        return (bound,pad) if bound>16*pad else None
    def handoff(start,end,query_wall,lower):
        proof=global_clearance(start,end,query_wall,lower)
        if proof is None:return wall_handoff_certificate(start,end,query_wall,lower)
        return True,dict(proof='GLOBAL_ORIGINAL_WALL_DISTANCE_MOTION_BOUND_MINUS_H_LOWER',
            minimum_g_nf_bound_m=proof[0],roundoff_m=proof[1],triangle_count=0)
    def geometric(start,end,query_wall,**kwargs):
        if global_clearance(start,end,query_wall,0.) is not None:
            return True,'GLOBAL_ORIGINAL_WALL_DISTANCE_SURFACE_MOTION_BOUND',None
        return swept_clearance_certificate(start,end,query_wall,**kwargs)
    interactions=bind_query_dependency(exact_interactions,{'wall_gap':cached_gap})
    assembly=bind_query_dependency(assemble_v1,{'exact_interactions':interactions})
    overrides={'wall_gap':cached_gap,'exact_interactions':interactions,'assemble_v1':assembly}
    if global_motion_bound:overrides.update(wall_handoff_certificate=handoff,swept_clearance_certificate=geometric)
    trial=bind_query_dependency(v1_trial,overrides)
    initial=bind_query_dependency(initial_world,{'wall_gap':cached_gap})
    method=bind_query_dependency(Particle65Stepper.step_to,{'initial_world':initial,'v1_trial':trial})
    return MethodType(method,stepper),gap_at.cache_info
