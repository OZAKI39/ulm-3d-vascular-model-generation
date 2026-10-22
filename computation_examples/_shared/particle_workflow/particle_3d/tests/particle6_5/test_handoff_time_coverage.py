import numpy as np

def test_no_skipped_or_duplicate_time(wall_run,pair_run,real_saved):
 for d in [wall_run,pair_run,real_saved]:
  rows=d['ledger'];assert rows[0]['t0_s']==0
  assert all(a['t1_s']==b['t0_s'] for a,b in zip(rows,rows[1:]))
  assert all(r['dt_s']>0 for r in rows)
  assert abs(sum(r['dt_s'] for r in rows)-d['summary']['horizon_s'])<=4*np.spacing(d['summary']['horizon_s'])
