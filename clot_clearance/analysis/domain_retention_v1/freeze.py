from common import *
from datetime import datetime,timezone
assert not (O/'DOMAIN_PROTOCOL_FREEZE.json').exists()
assert json.loads((O/'INITIAL_STATE_PARITY.json').read_text())['status']=='PASS'
assert json.loads((O/'FAILURE_SAFE_OUTPUT_AUDIT.json').read_text())['status']=='PASS'
old=json.loads((PREV/'EVENT_INTEGRATOR_FREEZE.json').read_text())
keys=['dt_base_s','minimum_h_s','checkpoint_predictor_trigger','damage_admission','material','gates','boundary','actual_analytic_expression','reference_sha256','ledger','normalization','budget','checkpoint','rollback','interval_strategy','precision','no_new_contact_calibration']
z={k:old[k] for k in keys}
for n,h in z['reference_sha256'].items():assert sha(O/'provenance'/n)==h
z.update(task=TASK,run_id=RUN,frozen_utc=datetime.now(timezone.utc).isoformat(),library_sha256=sha(LIB),driver_sha256={n:sha(A/n) for n in ['common.py','native.py','event_integrator.py','journal.py','run_protocol.py','preflight.py','freeze.py']},previous_freeze_sha256=sha(PREV/'EVENT_INTEGRATOR_FREEZE.json'),inherited_fields_exact=keys,domain=dict(bounds_m=[[-BOX,BOX]]*3,boundary='f f f',remap=False,no_wall=True),output_policy=dict(accepted_scalars_and_events='flush every accepted step; fsync at checkpoints and failure',node_stream='all accepted x/v/F in FP64',rolling_restart='two alternating slots every 40 accepted steps, plus initial/final and predicted clearance < dx',damage_restart='same previous trigger and rollback algorithm',failure='cached accepted FP64, no native query or restart after exception'),formal_cases=1,from_t0=True,end_time_s=.0007)
dump(O/'DOMAIN_PROTOCOL_FREEZE.json',z)
(O/'DOMAIN_PROTOCOL_FREEZE.sha256').write_text(sha(O/'DOMAIN_PROTOCOL_FREEZE.json')+'  DOMAIN_PROTOCOL_FREEZE.json\n')
print('FROZEN',sha(O/'DOMAIN_PROTOCOL_FREEZE.json'))
