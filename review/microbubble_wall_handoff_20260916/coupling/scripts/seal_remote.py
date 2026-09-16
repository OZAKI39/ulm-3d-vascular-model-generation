from pathlib import Path
import hashlib,json,shutil,subprocess,sqlite3
r=Path('/workspace/microbubble_lammps/results/palabos_lammps_coupling_v0_20260915_224019');w=Path('/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019')
assert json.loads((r/'EXECUTION_STATE.json').read_text())['status']=='PASS'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for f in (w/'provenance').glob('*'):
 if f.is_file():shutil.copy2(f,r/'provenance'/('PINNED_'+f.name))
for name in ['CMakeCache.txt','build.ninja']:
 p=w/'build'/name
 if p.exists():shutil.copy2(p,r/'provenance'/name)
con=sqlite3.connect(r/'cases/kokkos_compatibility/kokkos_host_fix_trace.sqlite');names=[n[0] for n in con.execute('select name from sqlite_master where type="table"') if 'MEMCPY' in n[0] and 'ENUM' in n[0]]
meaning={n:list(con.execute('select * from '+n)) for n in names};(r/'cases/kokkos_compatibility/CUPTI_ENUM_IDENTITIES.json').write_text(json.dumps(meaning,indent=2)+'\n');print(meaning)
(r/'provenance/BINARY_DYNAMIC_DEPENDENCIES.txt').write_text(subprocess.check_output(['ldd',str(w/'build/coupled_lmp_gpu')],text=True))
files=[p for p in sorted(r.rglob('*')) if p.is_file() and p.name not in ['SHA256SUMS','NUMERICAL_EVIDENCE_SHA256SUMS']]
text=''.join(sha(p)+'  '+str(p.relative_to(r))+'\n' for p in files)
(r/'NUMERICAL_EVIDENCE_SHA256SUMS').write_text(text)
(r/'SHA256SUMS').write_text(text+sha(r/'NUMERICAL_EVIDENCE_SHA256SUMS')+'  NUMERICAL_EVIDENCE_SHA256SUMS\n')
print('REMOTE_SEAL',sha(r/'SHA256SUMS'),'FILES',len(files)+1,'GiB',sum(p.stat().st_size for p in r.rglob('*') if p.is_file())/2**30)
