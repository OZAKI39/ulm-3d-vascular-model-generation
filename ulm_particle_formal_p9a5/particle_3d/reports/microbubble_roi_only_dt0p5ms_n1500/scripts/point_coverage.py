"""Descriptive point-flow diagnostic of the already frozen 1500 births.

These point trajectories do not enter the microbubble results or cohort
selection; finite-size paths remain the only reported microbubble trajectories.
"""
import json
import campaign

root,env,protected=campaign.setup()
from particle_3d.particle82_point_native import NativePointTracer
tracer=NativePointTracer(env)
events=json.loads((campaign.HERE/'data/cohort.json').read_text())['events']
rows=[]
for event in events:
    result=tracer.trace(event['birth_center_m'],step_m=.2e-6,error=1e-11,horizon_m=2e-3)
    rows.append(dict(particle_id=event['particle_id'],outlet=result['outlet'],
        end_reason=result['end_reason'],diameter_um=event['diameter_um'],
        role='POINT_FLOW_DIAGNOSTIC_ONLY_NOT_A_MICROBUBBLE_RESULT'))
campaign.dump(campaign.HERE/'data/point_birth_coverage.json',dict(rows=rows,
    cohort_fixed_before_diagnostic=True,used_for_selection=False,
    counts={o:sum(str(r['outlet'])==o for r in rows) for o in sorted({str(r['outlet']) for r in rows})}))
print('POINT_BIRTH_COVERAGE',json.dumps({o:sum(str(r['outlet'])==o for r in rows) for o in sorted({str(r['outlet']) for r in rows})}),flush=True)
