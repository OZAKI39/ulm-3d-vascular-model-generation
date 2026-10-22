"""Optional validation observer; default path installs nothing."""
from copy import deepcopy
from types import MethodType
import numpy as np
from .particle6_stepper import bind_query_dependency


def instrument(stepper, *, enabled=False, rows=None):
    if not enabled:return stepper
    if rows is None:raise ValueError('Explicit validation record sink required')
    method=stepper.advance_cached.__func__;factory=method.__globals__['v1_trial']
    original_solve=factory.__globals__['solve_resistance'];current=None
    def solve(system,**kwargs):
        result=original_solve(system,**kwargs)
        current.update(solver=deepcopy(result.record),free=system.free.tolist(),
            velocity=result.velocity.tolist(),hydro=result.unconstrained.tolist(),
            interactions=deepcopy(system.blocks),planar=deepcopy(getattr(system,'planar_diagnostics',[])))
        return result
    bound=bind_query_dependency(factory,{'solve_resistance':solve})
    def trial_factory(*args,**kwargs):
        trial=bound(*args,**kwargs)
        def observed(old,dt):
            nonlocal current
            current=dict(time_s=old.time_s,requested_dt_s=dt,
                position=next(iter(old.shapes.values())).center_m.tolist(),accepted=False)
            try:
                new,proof=trial(old,dt)
                current.update(accepted=True,accepted_dt_s=new.time_s-old.time_s,
                    handoff_event=new.projection.get('handoff_event'),
                    continuous_certificates=new.projection['continuous_certificates'])
                return new,proof
            except Exception as error:
                current['error']=str(error);raise
            finally:
                rows.append(current);current=None
        return observed
    stepper.advance_cached=MethodType(bind_query_dependency(method,{'v1_trial':trial_factory}),stepper)
    return stepper
