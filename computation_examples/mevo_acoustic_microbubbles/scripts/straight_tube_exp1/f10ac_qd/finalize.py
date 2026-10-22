"""Source/history preservation, artifact hashes and remote delivery identity."""
import time,json,subprocess,shutil
from pathlib import Path
from .common import *
def main():
 begin=time.perf_counter();verify();root=Path('/workspace/mevo_acoustic_microbubbles/f10ab_shape_20261018');delivery=read(root/'DELIVERY.json');source=Path(delivery['trusted_source_execution']);checks=[]
 for p in (ROOT/'src/mevo_mb').rglob('*.py'):
  rel=p.relative_to(ROOT);checks.append(dict(path=str(rel),matches=sha(p)==sha(source/rel)))
 history=[];manifest=read(AB/'INPUT_OUTPUT_SHA256.json')
 # Read the actual baseline manifest schema, preserving all available raw historical inputs.
 for p in AB.glob('*/**/trajectory.h5'):
  history.append(dict(path=str(p.relative_to(AB)),sha256=sha(p)))
 raw_manifest=read(AB/'RAW_DATA_MANIFEST.json')
 expected={r['raw_path']:r['raw_sha256'] for r in raw_manifest['records']}
 for r in history:r['matches_baseline']=r['path'] in expected and r['sha256']==expected[r['path']]
 assert len(history)==len(expected) and all(r['matches_baseline'] for r in history)
 write(BASE/'evidence/BASELINE_PRESERVATION.json',dict(delivery=delivery,source_checks=checks,all_source_matches=all(c['matches'] for c in checks),history=history,original_raw_manifest=raw_manifest,wall_s=time.perf_counter()-begin))
 assert all(c['matches'] for c in checks)
 shutil.copy2(ROOT/'scripts/straight_tube_exp1/f10ac_qd/README.md',BASE/'EXECUTION_AND_DATA.md')
 for p in Path('/workspace/mevo_acoustic_microbubbles/f10ac_ad_20261018/logs').iterdir():
  if p.is_file():shutil.copy2(p,BASE/'evidence/server_logs'/p.name)
 allfiles={}
 for base in [BASE,AD]:
  for p in base.rglob('*'):
   if p.is_file() and p.name!='INPUT_OUTPUT_SHA256.json' and not p.name.endswith('.lock'):allfiles[str(p.relative_to(ROOT))]=dict(sha256=sha(p),bytes=p.stat().st_size)
 write(BASE/'INPUT_OUTPUT_SHA256.json',dict(files=allfiles,identity=identity(),excluded=['self manifest','inactive locks']))
 server=Path('/workspace/mevo_acoustic_microbubbles/f10ac_ad_20261018');write(server/'DELIVERY.json',dict(identity=identity(),source_execution=str(ROOT),results=str(BASE),frequency_results=str(AD),manifest_sha256=sha(BASE/'INPUT_OUTPUT_SHA256.json'),training_ready=False))
if __name__=='__main__':main()
