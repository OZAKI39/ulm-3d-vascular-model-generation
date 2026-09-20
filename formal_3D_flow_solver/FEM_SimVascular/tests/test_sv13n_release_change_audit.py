from sv13n_support import *
def test_six_official_change_notes():
 d=accepted('release_change_audit');assert set(d['official_notes'])==set(map(str,range(320,326)))
 for n in d['official_notes'].values():assert hashlib.sha256((ROOT/n['path']).read_bytes()).hexdigest()==n['sha256']
 assert any('PetscBool' in n['symbols'] for n in d['relevant_notes'])
def test_old_and_new_real_runtime_view_formats():
 old=parse_runtime_semantics((R/'old_cpu_runtime_view_excerpt.txt').read_text())
 new=parse_runtime_semantics((ROOT/'logs/sv1_3n/remote/OFFICIAL_CPU_INITIAL.log').read_text())
 assert len(old)==1 and len(new)==4
 for d in old+new:runtime_semantics_gate(d)
