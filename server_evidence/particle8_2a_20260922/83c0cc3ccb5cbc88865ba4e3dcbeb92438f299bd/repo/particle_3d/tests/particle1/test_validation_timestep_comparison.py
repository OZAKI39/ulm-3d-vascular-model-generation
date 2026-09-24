import json
import numpy as np


def test_fixed_halving_and_saved_metrics_are_actual_trajectories(trajectory_cases,repo):
    summaries=[summary for _,summary in trajectory_cases]
    initialization=json.loads((repo/'particle_3d/reports/particle1/data/05_real_initialization.json').read_text())
    dt=[s['validation_dt_s'] for s in summaries]
    assert dt==initialization['validation_timesteps_s']
    assert dt[1]==dt[0]/2 and dt[2]==dt[0]/4
    for rows,summary in trajectory_cases:
        positions=np.array([[float(r[f'{a}_m']) for a in 'xyz'] for r in rows])
        measured=np.linalg.norm(np.diff(positions,axis=0),axis=1).sum()
        bound=32*np.finfo(float).eps*len(rows)*summary['trajectory_length_m']
        assert abs(measured-summary['trajectory_length_m'])<=bound
        assert float(rows[-1]['time_s'])==summary['exit_time_s']
        assert summary['timestep_role']=='VALIDATION_ONLY'
        assert summary['production_particle_timestep_frozen'] is False
