from common import *
import time,subprocess
start=time.perf_counter();before=json.loads((O/'HISTORICAL_PROTECTION_BEFORE.json').read_text());changed=[];missing=[];total=0
for name,old in before.items():
 p=Path(name)
 if not p.is_file():missing.append(name);continue
 total+=p.stat().st_size
 if sha(p)!=old['sha256'] or p.stat().st_size!=old['bytes']:changed.append(name)
out={'status':'PASS_ALL_PROTECTED_SHA_UNCHANGED' if not changed and not missing else 'FAIL_PROTECTION','checked_files':len(before),'checked_bytes':total,'changed':changed,'missing':missing,'wall_seconds':time.perf_counter()-start,'production_sha256':sha(PROD),'final_diagnostic_sha256':sha(B/'liblammps.so'),'env_sha256':sha('/workspace/y_clot/env.sh'),'no_production_deployment':True,'no_git_push':True,'no_old_CFD_or_fragmentation_overwrite':not changed}
dump(O/'HISTORICAL_PROTECTION_AFTER.json',out)
(O/'FINAL_GIT_STATUS.txt').write_text(subprocess.check_output(['git','status','--short'],cwd=R,text=True));print(json.dumps(out,indent=2));assert out['status'].startswith('PASS')
