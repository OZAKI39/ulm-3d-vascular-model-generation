"""Exact read-only query memoization; no interpolation or mechanics changes."""
from functools import lru_cache
from types import MethodType
import numpy as np

def install(env):
    from particle_3d.particle6_stepper import bind_query_dependency
    original_sample=env.field.sample
    @lru_cache(maxsize=32)
    def cached_sample(key):
        result=original_sample(np.frombuffer(key,dtype=np.float64))
        for value in vars(result).values():
            if isinstance(value,np.ndarray):assert not value.flags.writeable
        return result
    def sample(position):
        point=env.field._position(position)
        return cached_sample(point.tobytes())
    env.field.sample=sample
    # The entire wall is immutable. Repeating its maximum reduction does not
    # change the answer; preserve the original rounding expression verbatim.
    triangles=env.wall.triangles
    assert not triangles.flags.writeable
    maximum=float(np.max(np.abs(triangles)))
    original=env.wall.candidates.__func__
    roundoff=original.__globals__['roundoff_length']
    def cached_roundoff(*values):
        return roundoff(*(maximum if value is triangles else value for value in values))
    env.wall.candidates=MethodType(bind_query_dependency(original,{'roundoff_length':cached_roundoff}),env.wall)
    return cached_sample.cache_info
