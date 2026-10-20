from common import *
import difflib
def main():
 assert not (O/'LOAD_DIRECTION_FREEZE.json').exists()
 for name in ['PREPARED.json','NORMAL_MAPPING_AUDIT.json','LOAD_UNIT_CHECKS.json']:assert json.loads((O/name).read_text())['status']=='PASS'
 old=json.loads((B/'V3_PILOT_FREEZE.json').read_text());drivers={name:sha(A/name) for name in ['common.py','external.py','protocol_support.py','run_protocol.py']}
 frozen=dict(task=TASK,run_id=RUN,direction='REFERENCE_OUTWARD_CLOT_SURFACE_NORMAL',sign='OUTWARD_FROM_CLOT_INTO_LUMEN',scientific_label='SYNTHETIC_MECHANICAL_CONTROL_NOT_MICROBUBBLE_FORCE',historical_direction='HISTORICAL_C1_TANGENT',load_scalar_parameters=old['load'],same_input_force_magnitudes_equal=True,different_trajectories_can_change_scalar_forces=True,reference_normal_sha256=sha(O/'REFERENCE_NORMAL_FIELD.npz'),load_mapping_audit_sha256=sha(O/'NORMAL_MAPPING_AUDIT.json'),driver_sha256=drivers,historical_driver_sha256=old['driver_sha256'],input_sha256=old['input_sha256'],library_sha256=old['library_sha256'],gates=old['gates'],material_card=json.loads((B/'V3_MATERIAL_CARD.json').read_text()),boundary=old['boundary'],background_drag=old['background_drag'],energy=old['energy'],dt_base_s=BASE,minimum_substep_s=MIN_H,early_checkpoint_s=EARLY,end_s=END,physical_cases=1,checkpoint_continuation='Reload same accepted 0.8 ms native restart; preserve accumulated ledger; no accepted trajectory rerun',resource_policy=dict(MPI_ranks=1,threads=1,GPU=False,wall_limit=None,RSS_limit=None,historical_output_limit_GiB=2),candidate_method_sha256=sha(PREV_A/'graph_ops.py'),candidate_criteria=json.loads((PREV/'V3_CHUNK_CANDIDATES.json').read_text())['criteria'])
 dump(O/'LOAD_DIRECTION_FREEZE.json',frozen);(O/'LOAD_DIRECTION_FREEZE.sha256').write_text(sha(O/'LOAD_DIRECTION_FREEZE.json')+'\n')
 for name in ['external.py','run_protocol.py']:
  (O/'source_identity'/(name+'.diff')).write_text(''.join(difflib.unified_diff((HIST_A/name).read_text().splitlines(True),(A/name).read_text().splitlines(True),fromfile='historical/'+name,tofile='direction_control/'+name)))
 print('DIRECTION_FROZEN',sha(O/'LOAD_DIRECTION_FREEZE.json'),flush=True)
if __name__=='__main__':main()
