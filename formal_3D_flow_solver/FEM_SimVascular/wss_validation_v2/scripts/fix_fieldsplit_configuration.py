from pathlib import Path
import shutil,json,time,hashlib
V=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
c=V/'stage3/execution_checks/vessel_baseline_cpu8_fieldsplit';old=c.with_name(c.name+'_blocksize_conflict');assert not old.exists();assert 'Cannot change fieldsplit blocksize from 4 to 2' in (c/'run/solver.log').read_text();c.rename(old);c.mkdir();shutil.copytree(old/'SV_MESH',c/'SV_MESH');(c/'run').mkdir()
for name in ['solver.xml','PETSC_OPTIONS.txt']:shutil.copy2(old/'run'/name,c/'run'/name)
shutil.copy2(old/'policy.json',c/'policy.json');(c/'reports').mkdir();shutil.copy2(old/'reports/preregistration.json',c/'reports/preregistration.json')
s=(c/'run/PETSC_OPTIONS.txt').read_text();remove='-pc_fieldsplit_block_size 4 -pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 ';assert remove in s;(c/'run/PETSC_OPTIONS.txt').write_text(s.replace(remove,''))
files=[p for root in [c/'SV_MESH',c/'run'] for p in root.rglob('*') if p.is_file()]+[c/'policy.json'];put(c/'input_hashes.json',{str(p.relative_to(c)):sha(p) for p in files})
put(old/'reports/execution.json',dict(status='FAIL',no_qualified_result=True,reason='Experimental options duplicated native explicit velocity/pressure index sets and conflicted with native PC blocksize2; not a physical or WSS failure',solver_log_sha256=sha(old/'run/solver.log')))
put(c/'reports/configuration_fix.json',dict(unix=time.time(),source='petsc_impl.cpp747-780: native explicit IS u indices4*g+[0,1,2], p indices4*g+3; PCFieldSplitSetIS handles actual mapping',removed_options=remove,why='Use native index sets; remove redundant incompatible command-line strided fields. No solver source change.',previous_failed_case=str(old)))
print('Native velocity/pressure index sets now used without conflicting options.')
