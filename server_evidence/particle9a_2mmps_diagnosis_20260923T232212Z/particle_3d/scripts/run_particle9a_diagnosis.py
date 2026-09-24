#!/usr/bin/env python3
"""DIAGNOSTIC_ONLY bounded same-12 A/B and exact-prefix checkpoint replays."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
import gzip
import json
from pathlib import Path
import time
import numpy as np
from particle_3d.particle9a_diagnostics import Observer, wrap_nominal_steps
from particle_3d.particle9a_motion import Particle9AStepper
from particle_3d.particle81_simulation import SavedTrajectoryStepper, environment
from particle_3d.particle82a_integration import integrate_admitted
from particle_3d.particle6_stepper import bind_query_dependency
from particle_3d.particle9a_provenance import require_current_flow
from particle_3d.particle8_replay import canonical_hash


def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)): return [clean(v) for v in value]
    if isinstance(value, (np.bool_,)): return bool(value)
    if isinstance(value, (np.integer,)): return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


def dump(path, value):
    Path(path).write_text(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def run_job(job):
    event, model, output, checkpoint, logging = job
    if (Path(output)/'trajectories'/f'mb_{event["particle_id"]:06d}.json').exists():
        raise FileExistsError('Diagnostic evidence already exists; choose a new output directory')
    env = environment(); observer = Observer(env, event, model)
    created = []
    before = canonical_hash(event)
    def factory(*args, **kwargs):
        if model == 'P65':
            stepper = SavedTrajectoryStepper(*args, **kwargs)
        else:
            stepper = Particle9AStepper(*args, gradient_provider=lambda x:env.field.sample(x).velocity_gradient_s_inv, **kwargs)
        created.append(stepper)
        return wrap_nominal_steps(stepper, observer, enabled=logging, checkpoint_p65_time=checkpoint)
    integrate = bind_query_dependency(integrate_admitted, {'SavedTrajectoryStepper': factory})
    result = integrate(deepcopy(event), output=output)
    assert canonical_hash(event) == before
    result.update(diagnostic_only=True, model=model, logging_enabled=logging,
        dynamics='DIAGNOSTIC_ONLY_'+model, original_birth_metadata_sha256=before,
        flow=require_current_flow(), inlet_plane=observer.plane, checkpoint=observer.checkpoint,
        planar_wall_statistics=getattr(created[0], 'planar_statistics', None) if created else None,
        recorder_trial_count=len(observer.trials), recorder_state_count=len(observer.states))
    if checkpoint is not None and observer.checkpoint is None:
        raise AssertionError('Requested exact nominal checkpoint was not reached')
    folder = Path(output)/'trajectories'; stem = f'mb_{event["particle_id"]:06d}'
    dump(folder/(stem+'.json'), result)
    for kind, rows in [('trials', observer.trials), ('states', observer.states)]:
        with gzip.open(folder/(stem+'.'+kind+'.jsonl.gz'), 'wt') as f:
            for row in rows:
                f.write(json.dumps(clean(row), allow_nan=False, separators=(',', ':'))+'\n')
    return dict(particle_id=event['particle_id'], model=model, checkpoint=checkpoint,
        end_reason=result['end_reason'], failure_detail=result['failure_detail'],
        completed=result['completed'], seconds=result['wall_seconds'], samples=result['sample_count'],
        trials=len(observer.trials))


def main():
    p = argparse.ArgumentParser(); p.add_argument('--audit', required=True); p.add_argument('--output', required=True)
    p.add_argument('--workers', type=int, default=6); p.add_argument('--checkpoints')
    p.add_argument('--ids', type=int, nargs='+'); p.add_argument('--models', nargs='+', default=['P65', 'P9A'])
    p.add_argument('--diagnostics', action='store_true', help='DIAGNOSTIC_ONLY logging opt-in; default OFF')
    a = p.parse_args(); audit = json.loads(Path(a.audit).read_text())
    events = audit['events']; events = list(events.values()) if isinstance(events, dict) else events
    if a.ids: events = [e for e in events if e['particle_id'] in a.ids]
    assert len(events) <= 12 and len({e['particle_id'] for e in events}) == len(events)
    checkpoints = json.loads(Path(a.checkpoints).read_text()) if a.checkpoints else None
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    jobs = []
    for event in events:
        if checkpoints:
            if str(event['particle_id']) in checkpoints:
                if not a.diagnostics:raise ValueError('Checkpoint switch requires explicit --diagnostics')
                jobs.append((event, 'P9A_THEN_P65', str(out/'checkpoint'), float(checkpoints[str(event['particle_id'])]), a.diagnostics))
        else:
            for model in a.models:
                assert model in ('P65', 'P9A')
                jobs.append((event, model, str(out/model), None, a.diagnostics))
    dump(out/'run_config.json', dict(diagnostic_only=True, workers=a.workers, job_count=len(jobs),
        ids=[e['particle_id'] for e in events], event_hashes={e['particle_id']:canonical_hash(e) for e in events},
        flow=require_current_flow(), checkpoints=checkpoints, logging_enabled=a.diagnostics, production_invoked=False))
    start = time.monotonic(); results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(run_job, job) for job in jobs]):
            r = future.result(); results.append(r); print(json.dumps(r), flush=True)
            dump(out/'progress.json', results)
    dump(out/'completed.json', dict(results=results, wall_seconds=time.monotonic()-start))


if __name__ == '__main__': main()
