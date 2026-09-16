from pathlib import Path
import hashlib,json,subprocess
R=Path(__file__).resolve().parents[1];raw=R/'remote_raw';expected='5f6bdd8ecdb0527702e296ce12e80280246b060155c5995e83f301059618561d';actual=hashlib.sha256((raw/'SHA256SUMS').read_bytes()).hexdigest();assert actual==expected
names={line.split('  ',1)[1] for line in (raw/'SHA256SUMS').read_text().splitlines()};present={str(p.relative_to(raw)) for p in raw.rglob('*') if p.is_file() and p!=raw/'SHA256SUMS'};assert names==present,dict(extra=list(present-names),missing=list(names-present))
p=subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=raw,capture_output=True,text=True);(R/'logs/WSL_DOWNLOADED_SHA256_VERIFY.log').write_text(p.stdout+p.stderr);assert p.returncode==0
res=dict(status='PASS',manifest_sha256=actual,files_verified=len(names),exact_file_set_match=True,command='sha256sum -c SHA256SUMS',returncode=p.returncode);(R/'validation/REMOTE_TO_WSL_INTEGRITY.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
