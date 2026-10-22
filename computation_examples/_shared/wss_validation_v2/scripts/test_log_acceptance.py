from pathlib import Path
import json,sys,hashlib,re
V=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(V/'scripts/production'))
from log_acceptance import classify
base=V/'stage3/vessel_baseline_gpu_mpi1_halfdt';raw=(base/'run/solver.log').read_text();dt=json.loads((base/'policy.json').read_text())['dt_s']
tests=[('actual_recovered',raw,True),('missing_RHS_restore_proof',raw.replace('original_rhs_restored=1','original_rhs_restored=0'),False),('missing_zero_guess_proof',raw.replace('initial_guess_zero=1','initial_guess_zero=0'),False),('retry_not_healthy',re.sub(r'(SV13Q_END logical=3 attempt=1[^\n]*?)healthy=1',r'\g<1>healthy=0',raw),False),('premature_correction',raw.replace('SV13Q_BEGIN logical=3 attempt=1','NS 1-999 bogus correction\nSV13Q_BEGIN logical=3 attempt=1'),False)]
for name,expected in [('stage2/pipe_nr4',True),('stage3/vessel_baseline_gpu_mpi8',False),('stage3/execution_checks/vessel_baseline_cpu8_fieldsplit_blocksize_conflict',False)]:
 c=V/name;tests.append(('actual_'+c.name,(c/'run/solver.log').read_text(),expected))
rows=[]
for name,text,expected in tests:
 result=classify(text,dt);assert result['accepted']==expected,(name,result)
 rows.append(dict(test=name,expected_accept=expected,actual_accept=result['accepted'],failed_attempts=result['failed_KSP_attempts'],verified_recoveries=len(result['verified_recovered_attempts'])))
out=dict(all_tests_passed=True,tests=rows,classifier_sha256=hashlib.sha256((V/'scripts/log_acceptance.py').read_bytes()).hexdigest(),actual_recovered_log_sha256=hashlib.sha256((base/'run/solver.log').read_bytes()).hexdigest(),scope='Only validation monitor classification; raw failed attempt retained; no solver or tolerance change')
(V/'data/monitor_recovery_classification_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
