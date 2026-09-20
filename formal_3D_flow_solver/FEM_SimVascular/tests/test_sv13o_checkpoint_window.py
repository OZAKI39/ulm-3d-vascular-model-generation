from sv13o_support import *
def test_original_step60_checkpoint_and_native_reload():
 d=accepted('checkpoint_window');assert (d['start_step'],d['end_step'],d['timesteps'])==(60,70,10)
 cp=Path(d['checkpoint']['path']);assert hashlib.sha256(cp.read_bytes()).hexdigest()==d['runtime_checkpoint_sha256']
 assert d['checkpoint']['complete_integration_history'] and d['rank_count']==1 and not d['restart_layout_changed']
def test_all_windows_use_the_same_start_and_input():
 a=accepted('PERF_A_acceptance')
 for p in (R/'remote').glob('*_execution.json'):
  e=json.loads(p.read_text())
  if e['mode'] not in ('window','profile'):continue
  assert e['start_step']==60 and e['end_step']==70 and e['initial_checkpoint_sha256']==a['initial_checkpoint_sha256']
  assert e['config_sha256']==a['config_sha256'] and e['solver_sha256']==a['solver_sha256'] and e['PETSc_library_sha256']==a['PETSc_library_sha256']
  if result(e['name'])['status']=='PASS':assert {x['step'] for x in e['history']['linear_solves']}==set(range(61,71))

