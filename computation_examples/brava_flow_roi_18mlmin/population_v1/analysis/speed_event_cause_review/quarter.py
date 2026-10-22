"""The single pre-authorized additional dt control, for event822 only."""
import json
import signal
import numpy as np
from replay import HERE, OriginalRHS, save, run, exact_recorded_neighbor_check, threshold_bracket, plane_crossing


def main():
    previous=json.loads((HERE/'REPLAY_RESULTS.json').read_text())
    assert len(previous['runs'])==2
    reason=dict(event=822,decision='USE_ONE_RESERVED_QUARTER_DT_RUN',
        evidence='0.5 vs 0.25 ms: finite-interval increase differs by 1.957378 mm/s and full-window path by 422.006 micrometers',
        role='Resolve whether timestep differences materially change event attribution; not a convergence certification')
    save('QUARTER_DECISION.json',reason)
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('WINDOW_600S_LIMIT')))
    signal.alarm(600)
    ledger=json.loads((HERE/'ledger_822.json').read_text())
    with np.load(HERE/'original_822.npz') as z:original={k:z[k] for k in z.files}
    rhs=OriginalRHS();data,result=run(rhs,ledger,original,.000125)
    result['recorded_neighbors']=exact_recorded_neighbor_check(data)
    if not result['recorded_neighbors']['isolation_pass']:
        save('QUARTER_INVALID_NEIGHBOR.json',result)
        raise ValueError('ACTUAL_PAIR_INTERACTION_MISSING_FROM_QUARTER_REPLAY')
    common=original['physical_time_s'];ii=np.searchsorted(common,[.6505,.651]);at=data['common_indices']
    times=data['time_s'];speed=np.linalg.norm(data['RHS_m_s'],axis=1);local=(times>=.6505)&(times<=.653)
    event=dict(dt_upper_s=.000125,common_time_s=common[ii].tolist(),common_speed_mm_s=(speed[at[ii]]*1e3).tolist(),
        speed_increment_over_original_physical_interval_mm_s=float(np.diff(speed[at[ii]])[0]*1e3),
        maximum_adjacent_scalar_increment_mm_s=float(np.diff(speed)[local[:-1]&local[1:]].max()*1e3),
        first_10pct_bracket_s=threshold_bracket(times,speed,.3373128117608555+.1*.03622800676480964,.6505),
        first_90pct_bracket_s=threshold_bracket(times,speed,.3373128117608555+.9*.03622800676480964,.6505),
        original_event_plane_crossings=[plane_crossing(data,original['x_m'][k,1]) for k in ii])
    with np.load(HERE/'replay_822_dt250us.npz') as z:
        half_common_x=z['x_m'][z['common_indices']]
    result['maximum_path_delta_vs_half_at_common_m']=float(np.linalg.norm(data['x_m'][at]-half_common_x,axis=1).max())
    previous['runs'].append(result);previous['event_comparison'].append(event)
    previous['replay_count']=3
    previous['executed_engine_runs_including_initial_half_dt_without_saved_output']=4
    previous['additional_quarter_replay']=reason
    save('REPLAY_RESULTS.json',previous)
    print('QUARTER_COMPLETE',event,'path_vs_half',result['maximum_path_delta_vs_half_at_common_m'],flush=True)
    signal.alarm(0)


if __name__=='__main__':main()
