import json,hashlib,re

def test_review(report):
 text=(report/'PARTICLE7_REVIEW.md').read_text()
 assert len(re.findall(r'^## \d+\.',text,re.M))==13
 assert text.count('应该看什么：')==text.count('实际看到什么：')==text.count('有没有异常：')==15
 assert 'MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW' in text
 assert '静态障碍' in text and '10 nm' in text

def test_hashes_and_machine_results(report,repo):
 d=json.loads((report/'PARTICLE7_VALIDATION.json').read_text())
 for group,base in [('source_sha256',repo),('data_sha256',report),('figure_sha256',report)]:
  for name,digest in d[group].items(): assert hashlib.sha256((base/name).read_bytes()).hexdigest()==digest,name
 assert d['scheduled_mb_count']>=1000 and d['scheduled_rbc_count']>=2000
 assert d['admitted_mb_count']==d['scheduled_mb_count'] and d['admitted_rbc_count']==d['scheduled_rbc_count']
 assert d['restart_sequence_exact'] and d['duplicate_id_count']==d['lost_particle_count']==0
 assert d['birth_overlap_count']==d['birth_wall_violation_count']==d['birth_nearfield_handoff_violation_count']==0
 assert not d['particle8_started'] and d['no_cfd_executed']
 assert d['real_inlet_smoke']['accounting']['scheduled_mb_count']>0 and d['real_inlet_smoke']['accounting']['scheduled_rbc_count']>0
 assert d['manual_visual_review']=='PENDING_USER_REVIEW'

def test_long_independent_ledger_audit(report):
 d=json.loads((report/'data/11_independent_ledger_audit.json').read_text())
 assert d['status']=='PASS' and d['complete_csv_rows']>1000000
 assert d['all_birth_centers_on_inlet'] and d['mb_times_exact']
 assert d['wall_path_minimum_lower_bound_m']>0
