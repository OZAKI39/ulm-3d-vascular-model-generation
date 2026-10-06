"""Two bounded original-engine runs from ID822's actual recorded birth.

The other particles' original accepted event clock is held fixed. This is a
closed isolated-particle control, not full-population timestep convergence.
"""
from pathlib import Path
from types import SimpleNamespace
import json
import signal
import time
import numpy as np
from diagnose import HERE, OriginalRHS, save
from particle_3d.population_v1.gpu_engine import GPUPopulationEngine


class Recorder:
    def __init__(self):
        self.states = []
        self.audit = []
        self.terminals = []
    def snapshot(self, row):
        self.states.append(row)
    def _row(self, name, row):
        self.audit.append(row)
    def terminal(self, row):
        self.terminals.append(row)


def run(rhs, ledger, original, dt):
    initial = ledger['initial_state']
    assert initial['age'] == 0 and initial['birth_time'] == initial['physical_time_s']
    assert ledger['isolation_confirmed'] and not ledger['pair_rows']
    recorder = Recorder()
    engine = GPUPopulationEngine(rhs.env, dt=dt, pairs=True, realization_id=0,
        recorder=recorder, gpu_field=rhs.field, gpu_wall=rhs.wall)
    engine.time = initial['physical_time_s']
    particle = SimpleNamespace(particle_id=822, radius=initial['radius'],
        position=np.array(initial['x']), velocity=np.array(initial['v']),
        omega=np.array(initial['omega']), q=np.array(initial['q']))
    engine.insert_accepted(particle, ledger['birth'])
    t0, t1 = original['physical_time_s'][[0, -1]]
    ticks = np.arange(int(np.floor(t0/dt))+1, int(np.floor(t1/dt))+1)*dt
    targets = np.unique(np.r_[original['physical_time_s'][1:], ticks[(ticks>t0)&(ticks<t1)]])
    start = time.perf_counter()
    for target in targets:
        engine.advance(float(target))
        if not engine.pool.n:
            raise ValueError(('UNEXPECTED_REPLAY_TERMINAL', engine.terminals))
    engine.flush_output()
    rows = recorder.states
    times = np.array([r['physical_time_s'] for r in rows])
    x = np.array([r['x'][0] for r in rows])
    v = np.array([r['v'][0] for r in rows])
    assert np.all(np.diff(times)>0)
    at = np.searchsorted(times, original['physical_time_s'])
    if not np.array_equal(times[at], original['physical_time_s']):
        raise ValueError('MISSING_COMMON_ACCEPTED_PHYSICAL_TIMES')
    difference = np.linalg.norm(x[at]-original['x_m'], axis=1)
    # Between common endpoints both historical paths are linear; replay may
    # have extra substeps. Compare those points with the exact historical held
    # segment, not a newly integrated/interpolated CFD field.
    ix = np.searchsorted(original['physical_time_s'], times, side='right')-1
    ix = np.minimum(ix, len(original['physical_time_s'])-2)
    oldx = original['x_m'][ix] + (times-original['physical_time_s'][ix])[:,None]*original['velocity_m_s'][ix+1]
    bound = float(np.linalg.norm(x-oldx,axis=1).max())
    margin = ledger['minimum_swept_margin_above_pair_shell_m']
    sampled = rhs.batch(x, ledger['radius_m'])
    data = dict(time_s=times,x_m=x,held_m_s=v,RHS_m_s=sampled['U'][:,:3],
        tetra=sampled['tetra'],gap_m=sampled['wall'][:,0],actual_dt_s=np.diff(times),
        common_indices=at,path_delta_at_common_m=difference)
    np.savez_compressed(HERE/f'replay_822_dt{dt*1e6:.0f}us.npz',**data)
    result=dict(dt_upper_s=dt,elapsed_s=time.perf_counter()-start,
        actual_dt_s=dict(min=float(np.diff(times).min()),median=float(np.median(np.diff(times))),
            max=float(np.diff(times).max()),accepted_count=len(times)-1),
        maximum_path_delta_at_common_m=float(difference.max()),
        maximum_path_delta_at_all_replay_vertices_m=bound,
        remaining_neighbor_margin_m=margin-bound,
        backend_counts=dict(engine.backend_counts),audit=recorder.audit,terminals=recorder.terminals)
    save(f'replay_822_dt{dt*1e6:.0f}us.json',result)
    return data,result


def exact_recorded_neighbor_check(data):
    """Continuous swept pairs on the replay and other particles' stored paths.

    A failed global displacement bound is inconclusive. This exact check visits
    every other recorded particle, retains births, and splits at both clocks.
    """
    from diagnose import POP
    runroot=POP/'runs/realization_0000_contactfix'
    manifest=json.loads((runroot/'MANIFEST.json').read_text())
    times=data['time_s'];positions=data['x_m'];held=data['held_m_s']
    def replay_positions(q):
        idx=np.searchsorted(times,q,side='right')-1
        idx=np.clip(idx,0,len(times)-2)
        return positions[idx]+(q-times[idx])[:,None]*held[idx+1]
    previous=None;minimum=float('inf');closest=None;segments=0
    for chunk in manifest['chunks']:
        if chunk['last_time_s']<times[0] or chunk['first_time_s']>times[-1]:continue
        with np.load(runroot/chunk['path']) as z:
            ts=z['state_time_s'];offset=z['state_offsets'];idsall=z['particle_id'];xx=z['x'];rr=z['radius']
            for k,t in enumerate(ts):
                if t<times[0] or t>times[-1]:continue
                a,b=offset[k:k+2];ids=idsall[a:b];x=xx[a:b];r=rr[a:b]
                target=np.flatnonzero(ids==822)
                if not len(target):continue
                radius=r[int(target[0])]
                other=ids!=822
                pointgap=np.linalg.norm(x[other]-replay_positions(np.array([t]))[0],axis=1)-r[other]-radius
                shell=np.maximum(.05*np.minimum(r[other],radius),2e-9)
                value=float((pointgap-shell).min())
                if value<minimum:
                    minimum=value;closest=dict(other_id=int(ids[other][np.argmin(pointgap-shell)]),time_s=float(t),kind='ENDPOINT_INCLUDING_BIRTHS')
                if previous is not None:
                    oldt,oldids,oldx,oldr=previous
                    if t>oldt:
                        common,ia,ib=np.intersect1d(oldids,ids,return_indices=True)
                        ia,ib=ia[common!=822],ib[common!=822]
                        q=np.r_[oldt,times[(times>oldt)&(times<t)],t]
                        xp=replay_positions(q)
                        for j in range(len(q)-1):
                            alpha=(q[j:j+2]-oldt)/(t-oldt)
                            neighbors=oldx[ia][None,:,:]+alpha[:,None,None]*(x[ib]-oldx[ia])[None,:,:]
                            relative=neighbors-xp[j:j+2,None,:]
                            delta=relative[1]-relative[0];den=np.sum(delta*delta,axis=1)
                            frac=np.divide(-np.sum(relative[0]*delta,axis=1),den,out=np.zeros_like(den),where=den>0).clip(0,1)
                            shell=np.maximum(.05*np.minimum(r[ib],radius),2e-9)
                            gaps=np.linalg.norm(relative[0]+frac[:,None]*delta,axis=1)-r[ib]-radius-shell
                            value=float(gaps.min());segments+=len(ib)
                            if value<minimum:
                                minimum=value;l=int(np.argmin(gaps))
                                closest=dict(other_id=int(ids[ib[l]]),time_s=float(q[j]+frac[l]*(q[j+1]-q[j])),kind='CONTINUOUS_SWEPT_SEGMENT')
                previous=(float(t),ids.copy(),x.copy(),r.copy())
    return dict(minimum_margin_above_pair_shell_m=minimum,closest=closest,
                pair_segment_checks=segments,isolation_pass=minimum>0,
                role='All other historical held paths frozen; not a new full-population trajectory')


def threshold_bracket(t, speed, threshold, start):
    indices=np.flatnonzero((t>=start)&(speed>=threshold))
    if not len(indices):
        return None
    i=int(indices[0])
    return [float(t[max(0,i-1)]),float(t[i])]


def plane_crossing(data, y):
    # Exact geometry of the accepted held-velocity polyline, not reintegration.
    x,t=data['x_m'],data['time_s']
    ii=np.flatnonzero((x[:-1,1]-y)*(x[1:,1]-y)<=0)
    if not len(ii):return None
    i=int(ii[0]);fraction=(y-x[i,1])/(x[i+1,1]-x[i,1])
    return dict(time_s=float(t[i]+fraction*(t[i+1]-t[i])),
                accepted_bracket_s=[float(t[i]),float(t[i+1])])


def main():
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('WINDOW_600S_LIMIT')))
    signal.alarm(600)
    ledger=json.loads((HERE/'ledger_822.json').read_text())
    with np.load(HERE/'original_822.npz') as z:original={k:z[k] for k in z.files}
    rhs=OriginalRHS();runs=[];arrays=[]
    for dt in [.0005,.00025]:
        existing=HERE/f'replay_822_dt{dt*1e6:.0f}us.npz'
        if dt==.0005 and existing.exists():
            with np.load(existing) as z:data={k:z[k] for k in z.files}
            ds=np.diff(data['time_s'])
            result=dict(dt_upper_s=dt,status='REUSE_COMPLETED_NOMINAL_REPLAY_FROM_FIRST_ATTEMPT',
                evidence='monitoring/replay/job.log',elapsed_s=None,
                actual_dt_s=dict(min=float(ds.min()),median=float(np.median(ds)),max=float(ds.max()),accepted_count=len(ds)),
                maximum_path_delta_at_common_m=float(data['path_delta_at_common_m'].max()),
                backend_counts=None,backend_counts_reason='First attempt saved arrays and summary but failed the coarse bound on the subsequent half-dt run before its aggregate JSON; nominal integration not repeated')
        else:
            data,result=run(rhs,ledger,original,dt)
        result['recorded_neighbors']=exact_recorded_neighbor_check(data)
        if not result['recorded_neighbors']['isolation_pass']:
            save('REPLAY_INVALID_NEIGHBOR.json',result)
            raise ValueError('ACTUAL_PAIR_INTERACTION_MISSING_FROM_REPLAY')
        arrays.append(data);runs.append(result)
        print('CLOSED_REPLAY',dt,result['actual_dt_s'],'max_path_delta',result['maximum_path_delta_at_common_m'],flush=True)
    # The nominal rerun must reconstruct the actual event before attribution.
    if runs[0]['maximum_path_delta_at_common_m']>1e-10:
        raise ValueError('NOMINAL_REPLAY_DID_NOT_RECONSTRUCT_ORIGINAL_EVENT')
    common=original['physical_time_s'];ii=np.searchsorted(common,[.6505,.651])
    event=[]
    for data,result in zip(arrays,runs):
        at=data['common_indices'];v=data['RHS_m_s'];s=np.linalg.norm(v,axis=1)
        times=data['time_s'];s0,s1=s[at[ii]]
        local=(times>=.6505)&(times<=.653)
        event.append(dict(dt_upper_s=result['dt_upper_s'],common_time_s=common[ii].tolist(),
            common_speed_mm_s=(s[at[ii]]*1e3).tolist(),
            speed_increment_over_original_physical_interval_mm_s=float((s1-s0)*1e3),
            maximum_adjacent_scalar_increment_mm_s=float(np.diff(s)[local[:-1]&local[1:]].max()*1e3),
            first_10pct_bracket_s=threshold_bracket(times,s,.3373128117608555+.1*.03622800676480964,.6505),
            first_90pct_bracket_s=threshold_bracket(times,s,.3373128117608555+.9*.03622800676480964,.6505),
            original_event_plane_crossings=[plane_crossing(data,original['x_m'][k,1]) for k in ii]))
    save('REPLAY_RESULTS.json',dict(particle_id=822,status='CLOSED_ORIGINAL_GPU_ENGINE_REPLAY',
        role='DIAGNOSTIC_ONLY; other-particle historical event clock fixed',
        replay_count=2,executed_engine_runs_including_initial_half_dt_without_saved_output=3,
        runs=runs,event_comparison=event,
        actual_step_changed=runs[0]['actual_dt_s']['accepted_count']!=runs[1]['actual_dt_s']['accepted_count'],
        additional_quarter_replay='NOT_RUN: compare finite-interval signal with adjacent-sample amplitude; not full dt certification',
        other_targets='RHS reevaluation only; missing historical safety-monitor checkpoint state',
        limitation='No full-population re-integration or continuous-field physiological accuracy claim'))
    signal.alarm(0)


if __name__=='__main__':
    main()
