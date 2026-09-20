from sv13p_support import *
def test_required_decision_plots_exist_and_are_hashed():
 d=accepted('visuals');expected={'gpu_pc_candidates.png','gpu_pc_setup_cost.png','gpu_pc_apply_cost.png','gpu_pc_iterations.png','gpu_pc_memory.png'}
 if read('winner')['full_required']:expected|={'gpu_pc_winner_steady.png','gpu_pc_final_speedup.png'}
 assert {Path(f['path']).name for f in d['figures']}==expected
 for f in d['figures']:
  b=(ROOT/f['path']).read_bytes();assert b[:8]==b'\x89PNG\r\n\x1a\n' and len(b)>15000
  assert hashlib.sha256(b).hexdigest()==f['sha256']
