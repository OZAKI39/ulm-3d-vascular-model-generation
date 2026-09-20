#!/usr/bin/env python3
"""Full final old-project inventory, including .git, without optional Git writes."""
import json
import gzip
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,git_state,compare_inventory,write_json,sha256,now
R=ROOT/'reports/sv1';baseline=json.loads(gzip.decompress((R/'old_fem_baseline.json.gz').read_bytes()));old=Path(baseline['root'])
current=inventory(old);current['git']=git_state(old);(R/'old_fem_final.json.gz').write_bytes(gzip.compress((json.dumps(current,indent=2,ensure_ascii=False)+'\n').encode(),mtime=0))
result=compare_inventory(baseline,current);result.update(timestamp=now(),entries=len(current['files']),old_fem_root=str(old),
    git_state_unchanged=baseline['git']==current['git'],original_git=baseline['git'],current_git=current['git'],
    baseline_sha256=sha256(R/'old_fem_baseline.json.gz'),final_sha256=sha256(R/'old_fem_final.json.gz'),
    scope='Every old FEM regular file (including source, configs, tests, reports, outputs, logs, venv and .git) hashed, with mode/mtime/size and symlinks; atime excluded')
if not result['git_state_unchanged']:result['status']='FAIL'
write_json(R/'old_fem_readonly_audit.json',result)
print(json.dumps(result,indent=2));assert result['status']=='PASS'
