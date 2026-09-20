from sv13p_support import *
def test_native_checkpoint_science_and_window_frozen():
 policy=json.loads((ROOT/'configs/sv1_3p/policy.json').read_text())
 for d in cases():
  if d['mode']=='full':continue
  assert d['initial_checkpoint_sha256']==policy['checkpoint_sha256'] and d['start_step']==60
  path=ROOT/'outputs/sv1_3p'/d['name']/'solver.xml'
  scientific_xml_gate(ROOT/'outputs/sv1_3o/REAL_VASCULAR_GPU_PERF/solver.xml',path)
  if d['status']=='PASS' and d['mode']=='window':assert d['steps']==10 and d['stop_step']==70
 assert not list((ROOT/'configs/sv1_3p/runplans').glob('*BASELINE*WINDOW*'))
