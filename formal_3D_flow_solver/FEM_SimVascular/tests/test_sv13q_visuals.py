from sv13q_support import *
def test_required_decision_figures_exist():
 expected={'ilu_rebuild_candidates.png','ilu_rebuild_count.png','ilu_iterations_tradeoff.png','ilu_setup_time.png','early_vs_late_window.png','adaptive_rebuild_timeline.png'}
 if read('winner')['full_required']:expected|={'ilu_winner_steady.png','ilu_final_speedup.png'}
 d=read('visuals');assert d['status']=='PASS';assert set(d['files'])==expected
 for n in expected:
  p=R/n;assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n') and p.stat().st_size>5000
