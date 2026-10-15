"""Run in WSL only, after scientific PASS. Preserve complete previous bytes."""
from pathlib import Path
import hashlib,json,difflib,os,shutil
from datetime import datetime,timezone
root=Path('/home/lzy/projects/temp_storage/clot_github_sync_20261011/clot_clearance')
out=root/'verification/lps_native_softening_domain_retention_and_failure_safe_output_v1/20261014T150825Z'
p=Path('/home/lzy/projects/clot_rupture.md');folder=out/'wsl_document';folder.mkdir(exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
assert not (out/'WSL_DOCUMENT_UPDATE.json').exists(),'Document update already recorded'
before=p.read_bytes();text=before.decode();insert=(out/'WSL_DOCUMENT_INSERT.md').read_text();marker='# ⑥ 外加剪切'
assert '本轮计算域保留与异常保存审查（20261014T150825Z' not in text
assert text.count(marker)==1
for word in ['0.07170070833333','0.06969623734234','173.156623','230.776817',*list('①②③④⑤⑥⑦⑧⑨')]:assert word in text
index=text.index(marker);after=(text[:index]+insert+text[index:]).encode()
assert after.decode().replace(insert,'',1)==text
(folder/'clot_rupture.before.md').write_bytes(before)
(folder/'clot_rupture.after.md').write_bytes(after)
(folder/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),after.decode().splitlines(True),fromfile=str(p)+'.before',tofile=str(p)+'.after')))
assert sha(p.read_bytes())==sha(before),'Concurrent document edit: stop before write'
tmp=p.with_name(p.name+'.domain_retention.tmp');shutil.copy2(p,tmp)
with tmp.open('wb') as f:f.write(after);f.flush();os.fsync(f.fileno())
os.replace(tmp,p);actual=p.read_bytes();assert actual==after
record=dict(status='PASS',actual_WSL_path=str(p),utc=datetime.now(timezone.utc).isoformat(),before_sha256=sha(before),after_sha256=sha(after),actual_readback_sha256=sha(actual),readback_exact=True,all_historical_text_preserved_by_insertion=True,chapters_1_to_9_preserved=True,G0_nominal_values_preserved=True,before_bytes=len(before),after_bytes=len(after),insert_start_line=text[:index].count('\n')+1,insert_lines=insert.count('\n'),backup='wsl_document/clot_rupture.before.md',after='wsl_document/clot_rupture.after.md',diff='wsl_document/clot_rupture.diff')
(out/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n');(folder/'ACTUAL_WSL_READBACK.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(record,ensure_ascii=False,indent=2))
