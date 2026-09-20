from sv13o_support import *
def test_five_development_plots_are_real_and_hashed():
 d=accepted('visuals');assert {Path(f['path']).name for f in d['figures']}=={'gpu_tuning_candidates.png','ksp_cost_comparison.png','gpu_time_breakdown.png','optimized_steady_convergence.png','gpu_runtime_before_after.png'}
 for f in d['figures']:
  p=ROOT/f['path'];assert p.stat().st_size>10000 and hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256']
