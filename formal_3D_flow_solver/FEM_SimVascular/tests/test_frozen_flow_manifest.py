import math
from fem_sync_support import ROOT, read, digest


def test_frozen_flow_is_actual_stage_q_final_state():
    f=read('frozen_reference/flow/flow_field_manifest.json')
    a=read('reports/sv1_3q/REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance.json')
    assert digest(ROOT/f['path'])==f['sha256']==a['reload']['sha256']
    assert f['step']==a['stop_step']
    assert f['velocity']==dict(name='Velocity',association='POINT',components=3,dtype='float64',units='m/s',finite=True)
    assert f['pressure']['name']=='Pressure' and f['pressure']['association']=='POINT'
    assert f['cell_arrays']==[]
    for name in ['velocity_gradient','vorticity','strain_rate']:
        assert f[name]=='NOT STORED IN FROZEN VTU'


def test_same_run_time_config_checkpoint_and_safe_stop():
    f=read('frozen_reference/flow/flow_field_manifest.json')
    a=read('reports/sv1_3q/REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance.json')
    s=read('reports/sv1_3q/steady_history.json')
    c=read('frozen_reference/fem_checkpoint/checkpoint_manifest.json')
    assert math.isclose(f['TimeValue'],c['time_s'],rel_tol=1e-12)
    assert f['safe_stop_step']==c['step']==s['stop_step']
    assert f['monitor_last_saved_step']==f['first_formal_steady_step']==s['first_full_steady_step']
    assert digest(ROOT/f['solver_config']['path'])==a['config_sha256']
    assert digest(ROOT/c['path'])==c['sha256']


def test_preserved_performance_counts_are_observed():
    b=read('frozen_reference/baseline_summary.json')
    a=read('reports/sv1_3q/REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance.json')
    assert b['statistics']==a['statistics']
    assert b['ILU_rebuild_count']==a['profile']['events']['MatLUFactorNum']['count']
    assert b['ILU_rebuild_count']+b['ILU_reuse_count']==b['statistics']['KSP_solves']
    assert b['linear_unrecovered_failures']==b['nonlinear_failures']==0
