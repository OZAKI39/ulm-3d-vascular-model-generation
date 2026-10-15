"""Immutable preregistration, before new F10S physics evaluations."""
from dataclasses import replace
from datetime import datetime, timezone
import numpy as np
from mevo_mb.control.scheduled_pzt import Schedule
from .common import *

def main():
    assert not (BASE/'PREREGISTRATION.json').exists()
    p=dict(baseline_commit='9f1db411836e9a67dff1f77a344e5457674e6a81',labels=read(R/'PREREGISTRATION.json')['protocol']['labels'],
       training_ready=False, hardware_executed=False, frequency_Hz=2.5e6,flow='STRICT_ZERO',
       total_scientific_budget_s=600,remote_budget_s=380,per_run_budget_s=65,max_trajectories=11,
       route_budget_s={'A':200,'B':180,'C':100,'D':120},
       target_m=[2e-7,0,0],target_tolerance_m=7.5e-8,strict_endpoint_tolerance_m=1e-8,
       target_definition='Absolute world coordinate; fixed before runs',A_duration_s=.003,B_screen_s=.0005,B_duration_s=.002,
       B_source_Pa_equivalent=400,source_limit_Pa_equivalent=1500,dt_s=1e-6,dt_half_s=.5e-6,cache_spacing_m=[3e-6,1.5e-6],
       initial='Original fixed-ID 32-particle manual contact demo. Every run freshly initialized. Not acoustic formation.',
       admission_policy=policy(), C=dict(amplitudes_Pa_equivalent=[400,1036.2016581041792],phases_deg=[-6,0,6],
       peak_search_halfwidth_m=70e-6,force_search_halfwidth_m=8e-6,jacobian_delta_m=50e-9,
       phase_step_deg=1,frozen_diagnostic_scope='FROZEN_SHAPE_COM_DIAGNOSTIC; fixed internal geometry; not physical rigid-body validation',
       long_time_stop='No full moving trap if relaxation time >2ms or compression risk. No long trajectories.'),
       D=dict(Y_distance_offsets_m=[-25e-6,-5e-6,0,5e-6,25e-6],COM_y_offsets_m=[-5e-6,0,5e-6],
       frequencies_Hz=[2.475e6,2.5e6,2.525e6],single_phase_deg=[0,90,180],source_Pa_equivalent=1500,
       scope='HOMOGENEOUS_WATER_PROXY_SPECIFIC; exploratory geometry offsets, not measured assembly error',
       off_frequency='FIELD_ONLY; no APECSS forces or trajectories outside 2.5MHz',
       pdms_preflight=dict(PPW=[6,8],CFL=.2,PML_cells=20,water_envelope_m=.001,gpu_budget_GiB=20,
       frequency_Hz=2.5e6,material='Xu2020 SYLGARD184 reference measured 3-7MHz; 2.5MHz not validated; no silent extrapolation')))
    write(CONFIG,p)
    old=Schedule.from_record(read(R/'navigation_v1/runs/X_PAIR_DT1/INPUT.json')['schedule'])
    pair=Schedule((replace(old.segments[0],duration_s=.003),))
    seq=Schedule.from_record(read(R/'RUN_SCHEDULES.json')['ONE_HOT_X_SEQUENCE'])
    plans={'X_PAIR_3MS':pair.record()}
    for key,segment in zip(['X_MINUS_ONLY','X_PLUS_ONLY'],seq.segments):
        assert max(segment.amplitude_Pa_equivalent)==400
        plans[key]=Schedule((replace(segment,segment_id=0,start_s=0,duration_s=.002),)).record()
    write(BASE/'RUN_SCHEDULES.json',plans)
    np.savez_compressed(BASE/'INITIAL_STATES.npz',original=np.load(R/'INITIAL_STATES.npz')['original'])
    files=[CONFIG,BASE/'RUN_SCHEDULES.json',BASE/'INITIAL_STATES.npz',R/'navigation_v1/runs/X_PAIR_DT1/INPUT.json',R/'navigation_v1/runs/X_PAIR_DT1/trajectory.h5']
    files += [OLD/n for n in read(Q/'PREREGISTRATION.json')['input_sha256']]
    write(BASE/'PREREGISTRATION.json',dict(registered_utc=datetime.now(timezone.utc).isoformat(),protocol=p,input_sha256={str(f.relative_to(ROOT)):sha(f) for f in files},identity=identity()))
    print('F10S A/B/C/D registered; full-precision F10R controls and original IDs retained')

if __name__=='__main__':main()
