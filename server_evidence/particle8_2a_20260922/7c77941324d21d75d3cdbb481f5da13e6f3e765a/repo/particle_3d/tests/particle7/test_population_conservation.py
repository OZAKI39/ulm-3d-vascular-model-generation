import pytest
def test_conservation(engine):
 row=engine.step_to(.005)
 for species in ['mb','rbc']:
  assert row['scheduled_'+species+'_count']==row['admitted_'+species+'_count']+row['pending_'+species+'_count']
  assert row['admitted_'+species+'_count']==row['active_'+species+'_count']+row['exited_'+species+'_count']
 assert row['scheduled_rbc_volume_m3']==pytest.approx(row['active_rbc_volume_m3']+row['exited_rbc_volume_m3'],rel=1e-14)

def test_saved_long_prefix_replays_current_code(report,tmp_path):
 import csv,gzip,itertools
 from particle_3d.particle7_long import run_long
 run_long(tmp_path,mb_target=1)
 with gzip.open(tmp_path/'11_long_lifecycle.csv.gz','rt') as short,gzip.open(report/'data/11_long_lifecycle.csv.gz','rt') as saved:
  a=list(itertools.islice(csv.DictReader(short),1200)); b=list(itertools.islice(csv.DictReader(saved),1200))
 assert len(a)==len(b)==1200 and a==b
 assert any(e['species']=='MB' for e in a)

def test_tube_hct_is_uniform_time_diagnostic(report):
 import json,numpy as np
 d=json.loads((report/'data/14_tube_hct_statistics.json').read_text())
 rows=np.genfromtxt(report/'data/14_uniform_time_tube_hct.csv',delimiter=',',names=True)
 assert np.ptp(np.diff(rows['time_s']))<1e-12
 assert np.std(rows['observed_tube_hct'])>.001 and rows['observed_tube_hct'][0]==0
 assert d['sampling']=='UNIFORM_PHYSICAL_TIME_NOT_BIRTH_CONDITIONED' and not d['enforced']
 # With equal plug velocity and no backlog, residence-time conservation predicts
 # the feed fraction apart from the explicitly empty initial control section.
 assert abs(d['exact_residence_time_weighted_mean']-.45)<1e-5
 assert d['real_lumen_snapshot']!=.45
