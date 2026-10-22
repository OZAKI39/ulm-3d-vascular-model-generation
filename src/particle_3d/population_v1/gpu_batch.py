"""Batch independent realizations without changing any realization's time grid.

Only the numerical GPU work is stacked. Each realization retains its own exact
birth/event/refinement clock. A neighbor event in realization 2 never inserts
an extra Euler subdivision in realization 1 (which would change its solution).
"""
import numpy as np
import torch
from .gpu_resistance import assemble_isolated,solve_isolated


def prefetch(engines):
    live=[e for e in engines if e.pool.n]
    if not live:return dict(realizations=0,particles=0)
    if len({e.realization_id for e in live})!=len(live):raise ValueError('Distinct realization IDs required')
    if any(e.env is not live[0].env for e in live):raise ValueError('Batch must share one authoritative flow environment')
    x=torch.cat([e.pool.x[:e.pool.n] for e in live]);radius=torch.cat([e.pool.radius[:e.pool.n] for e in live])
    hint=torch.cat([e.pool.tetra_id[:e.pool.n] for e in live]);field=live[0].field;wall=live[0].wall
    sampled=field.sample(x,hint)
    if bool((sampled[0]<0).any()):raise ValueError('BATCH_CENTER_OUTSIDE_FEM')
    if wall is not None:geometry=wall.query(x,radius)
    else:
        w=torch.zeros((len(x),14),device=x.device,dtype=x.dtype);w[:,0]=1.;w[:,7]=1.
        geometry=(torch.full((len(x),),-1,device=x.device,dtype=torch.int64),w,None,
            torch.zeros((len(x),32,3),device=x.device,dtype=x.dtype),torch.zeros(len(x),device=x.device,dtype=torch.int64))
    data=sampled[1];free=torch.cat([data[:,:3],.5*data[:,13:16]],dim=1)
    matrix,rhs,diagonal,assembly=assemble_isolated(radius,live[0].env.mu,geometry[1],data[:,4:13].reshape(-1,3,3),free)
    velocity,hydro,audit=solve_isolated(matrix,rhs,diagonal,free,geometry[3],geometry[4])
    offset=0
    for e in live:
        sl=slice(offset,offset+e.pool.n);offset+=e.pool.n
        cut=lambda value:value[sl] if isinstance(value,torch.Tensor) else value
        e.prefetched=dict(time=e.time,count=e.pool.n,sampled=tuple(cut(a) for a in sampled),
            wall=tuple(cut(a) for a in geometry),velocity=velocity[sl].clone(),hydro=hydro[sl],
            audit={k:cut(v) for k,v in audit.items()},assembly={k:cut(v) for k,v in assembly.items()})
    return dict(realizations=len(live),particles=len(x),single_batched_fem_wall_and_resistance_dispatch=True)


def perfuse_batch_until(engines,sources,end_times):
    if not (len(engines)==len(sources)==len(end_times)):raise ValueError('One source/end time per realization required')
    dispatches=[]
    while True:
        work=[]
        for e,s,end in zip(engines,sources,end_times):
            if e.time>=end:continue
            row=s.peek();ratio=e.time/e.dt
            tick=(int(np.floor(ratio+8*abs(np.spacing(ratio))))+1)*e.dt
            target=min(float(end),tick,row['birth_time_s'] if row is not None else float('inf'))
            work.append((e,s,target,row))
        if not work:break
        dispatches.append(prefetch([w[0] for w in work]))
        for e,s,target,row in work:
            e.advance(target)
            if row is not None and e.time==row['birth_time_s']:s.admit(e)
    return dict(dispatches=dispatches,accounting=[e.accounting() for e in engines])
