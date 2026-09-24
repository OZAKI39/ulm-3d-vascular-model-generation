"""Zero-radius P1 advection; no MB wall force, radius or handoff.

Native double-precision tetra-adjacency RK23 supplies the current diagnostic.
The initial VTK RK45 attempt remains available as a documented failed control.
All unresolved seeds stay in the denominator. Original cap triangles audit
actual local ODE exit segments; no outlet is assigned by nearest distance.
"""
from pathlib import Path
import json
import hashlib
import multiprocessing as mp
import os
import socket
import time
import numpy as np
from .particle82_provenance import atomic_json, sha256, require_remote
from .particle8_replay import canonical_hash

_ENV = None


def init_environment():
    global _ENV
    from .particle81_simulation import environment
    from .validation_boundary import ValidationBoundaryClassifier
    _ENV = environment()
    _ENV.point_boundaries = ValidationBoundaryClassifier({k:v for k,v in _ENV.boundaries.items() if k != 'INLET'})
    return _ENV


def trace_positions(positions, *, step_m=0.2e-6, error=1e-11, horizon_m=2e-3):
    from .particle82_point_native import NativePointTracer
    env = _ENV or init_environment()
    if not hasattr(env,'native_point'):env.native_point=NativePointTracer(env)
    rows=[]
    for i,p in enumerate(positions):
        row=env.native_point.trace(p,step_m=step_m,error=error,horizon_m=horizon_m)
        row['seed_index']=i;rows.append(row)
    return rows


def trace_positions_vtk(positions, *, step_m=0.2e-6, error=1e-9, horizon_m=2e-3):
    import pyvista as pv
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    env = _ENV or init_environment()
    mesh = env.field._grid.copy(deep=False)
    mesh.point_data['Velocity'] = env.field.velocity_nodes_m_s
    mesh.set_active_vectors('Velocity')
    seeds = pv.PolyData(np.asarray(positions, dtype=np.float64))
    tracer = vtk.vtkStreamTracer()
    tracer.SetInputData(mesh); tracer.SetSourceData(seeds)
    tracer.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, 'Velocity')
    tracer.SetIntegrationDirectionToForward(); tracer.SetIntegratorTypeToRungeKutta45()
    tracer.SetIntegrationStepUnit(vtk.vtkStreamTracer.LENGTH_UNIT)
    tracer.SetInitialIntegrationStep(step_m); tracer.SetMaximumIntegrationStep(step_m)
    tracer.SetMinimumIntegrationStep(step_m/128); tracer.SetMaximumError(error)
    tracer.SetMaximumPropagation(horizon_m); tracer.SetMaximumNumberOfSteps(40000)
    tracer.SetTerminalSpeed(1e-15); tracer.SetComputeVorticity(False)
    tracer.Update()
    output = tracer.GetOutput()
    rows = [dict(seed_index=i, outlet=None, end_reason='VTK_NO_PATH', vtk_reason=None,
                 path=np.array([[0., *p]], dtype=np.float64)) for i,p in enumerate(positions)]
    if output.GetNumberOfCells() == 0: return rows
    points = vtk_to_numpy(output.GetPoints().GetData())
    times = vtk_to_numpy(output.GetPointData().GetArray('IntegrationTime'))
    seed_ids = vtk_to_numpy(output.GetCellData().GetArray('SeedIds'))
    reasons = vtk_to_numpy(output.GetCellData().GetArray('ReasonForTermination'))
    for ci in range(output.GetNumberOfCells()):
        cell = output.GetCell(ci); ids = np.array([cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())])
        seed = int(seed_ids[ci]); path = np.column_stack([times[ids], points[ids]])
        reason = 'UNCLASSIFIED_VTK_TERMINATION'; outlet = None; closure = []
        # The last RK45 interior point can be less than one spatial step from
        # an open cap. Continue the ODE locally; never assign an outlet by distance.
        for attempt in range(4):
            p = path[-1,1:]; sample = env.field.sample(p)
            if not sample.inside_lumen:
                reason = 'LOCATOR_OUTSIDE_AT_SAVED_POINT'; break
            speed = np.linalg.norm(sample.velocity_m_s)
            if speed <= 1e-15: reason = 'ZERO_SPEED'; break
            dt = step_m/speed
            end = p + dt*sample.velocity_m_s
            hit = env.point_boundaries.first_event(p, end)
            if hit is not None:
                end = hit.position_m; dt *= hit.segment_fraction
                if dt > 0: path = np.vstack([path, [path[-1,0]+dt,*end]])
                closure.append(dict(dt_s=float(dt),velocity_m_s=sample.velocity_m_s.tolist(),tetra_id=sample.tetra_id))
                outlet = hit.role if hit.role.startswith('OUTLET_') else None
                reason = hit.role if outlet else 'WALL_CROSSING_NUMERICAL_DIAGNOSTIC_FAILURE'
                break
            if not env.field.sample(end).inside_lumen:
                reason = 'UNCLASSIFIED_DOMAIN_EXIT'; break
            path = np.vstack([path,[path[-1,0]+dt,*end]])
            closure.append(dict(dt_s=float(dt),velocity_m_s=sample.velocity_m_s.tolist(),tetra_id=sample.tetra_id))
        rows[seed] = dict(seed_index=seed,outlet=outlet,end_reason=reason,vtk_reason=int(reasons[ci]),
                          path=path,terminal_ode_steps=closure)
    return rows


def shard_job(task):
    first, positions, folder, config, provenance = task
    require_remote(provenance, provenance['hostname'])
    folder = Path(folder); stem=f'tracer_{first:07d}_{first+len(positions)-1:07d}'
    meta_path=folder/(stem+'.json'); path_file=folder/(stem+'.npz')
    seed_sha=hashlib.sha256(np.asarray(positions,dtype=np.float64).tobytes()).hexdigest()
    if meta_path.exists():
        meta=json.loads(meta_path.read_text())
        if meta['path_sha256'] != sha256(path_file): raise ValueError('Corrupt point-tracer shard')
        if meta.get('config_sha256')!=canonical_hash(config) or meta.get('seed_positions_sha256')!=seed_sha or meta.get('git_commit')!=provenance['source_git_commit']:
            raise ValueError('Point-tracer resume source/config/seed mismatch')
        if meta['hostname']!=provenance['hostname']:raise ValueError('Point-tracer shard host mismatch')
        return meta
    start=time.time(); rows=trace_positions(positions,**config)
    offsets=np.cumsum([0]+[len(r['path']) for r in rows])
    packed=np.concatenate([r.pop('path') for r in rows])
    temp=path_file.with_suffix('.tmp.npz')
    np.savez_compressed(temp,positions=positions,offsets=offsets,paths=packed)
    os.replace(temp,path_file)
    for r in rows:r['tracer_id']=first+r.pop('seed_index')
    meta=dict(role='ZERO_RADIUS_POINT_TRACERS_DIAGNOSTIC_ONLY',first_id=first,count=len(rows),rows=rows,
        hostname=socket.gethostname(),pid=os.getpid(),start_time=start,end_time=time.time(),
        worker_id=mp.current_process().name,git_commit=provenance['source_git_commit'],
        seed_positions_sha256=seed_sha,
        frozen_input_sha256=provenance['frozen_input_sha256'],config_sha256=canonical_hash(config),
        config=config,path_sha256=sha256(path_file),path_file=path_file.name)
    atomic_json(meta_path,meta)
    return meta


def run(positions, output, provenance, workers=8, shard_size=512, config=None):
    require_remote(provenance,provenance['hostname'])
    folder=Path(output);folder.mkdir(parents=True,exist_ok=True)
    init_environment()  # immutable mesh/field inherited by Linux fork (COW)
    config=config or dict(step_m=0.2e-6,error=1e-11,horizon_m=2e-3)
    from .particle82_point_native import NativePointTracer
    _ENV.native_point=NativePointTracer(_ENV)
    tasks=[(i,np.asarray(positions[i:i+shard_size]),str(folder),config,provenance) for i in range(0,len(positions),shard_size)]
    result=[]
    with mp.get_context('fork').Pool(workers) as pool:
        for meta in pool.imap_unordered(shard_job,tasks,chunksize=1):
            result.append(meta)
            print('POINT_TRACERS',sum(r['count'] for r in result),'/',len(positions),flush=True)
    return sorted(result,key=lambda r:r['first_id'])
