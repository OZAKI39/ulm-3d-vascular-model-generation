"""Explicit configured bindings around the unchanged finite-size P9A1/P6.5 integrator."""
from functools import partial
from . import formal_dynamics_p9a5 as original
from .formal_cohort_p9a5 import HorizonSchedule
from .particle6_stepper import bind_query_dependency
from .particle9a_motion import Particle9AStepper

def bind_job(env, identity, dt, *, accelerated=False):
    """Keep dt, force/contact algorithms and 3/6/12 s horizon ordering from the chosen case."""
    audit=bind_query_dependency(original.audit_result,{'ENV':env})
    job=bind_query_dependency(original.job,dict(ENV=env,IDENTITY=identity,DT=dt,
        audit_result=audit,HorizonSchedule=partial(HorizonSchedule,dt=dt)))
    if accelerated:
        from .native_adapter import accelerated_job,native_closest
        native_closest(env.wall.triangles[0,0],env.wall.triangles[:1])
        job=accelerated_job(job)
    return job
