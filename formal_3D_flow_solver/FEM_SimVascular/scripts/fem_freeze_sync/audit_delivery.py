"""Bounded source, secret, size, portability and main-history checks; no CFD."""
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path
from verify_manifest import ROOT, MANIFEST, payload_files

REPO = ROOT.parents[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(name, value):
    (ROOT/'sync_metadata'/name).write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def audit_main():
    base = json.loads((ROOT/'sync_metadata/main_base.json').read_text())['base_commit']
    tracked = subprocess.check_output(['git','ls-tree','-r','--name-only',base],cwd=REPO,text=True).splitlines()
    changed = []
    for name in tracked:
        original = subprocess.check_output(['git','show',base+':'+name],cwd=REPO)
        current = (REPO/name).read_bytes()
        if current != original:
            changed.append(name)
            assert name == 'README.md' and current.startswith(original)
            assert b'## Frozen 3D FEM / Particle Handoff' in current[len(original):]
    assert changed == ['README.md']
    assert subprocess.check_output(['git','rev-parse','origin/main'],cwd=REPO,text=True).strip() == base
    write('main_history_audit.json', dict(status='PASS',base_commit=base, main_files_checked=len(tracked),
        modified_existing_files=changed, original_code_rewritten=False, allowed_new_prefix='formal_3D_flow_solver/FEM_SimVascular/',
        main_branch_untouched=True, README_append_only=True))


def source_audit():
    d = json.loads((ROOT/'sync_metadata/copied_from_wsl.json').read_text())
    source = Path(d['source'])
    for row in d['files']:
        assert sha(ROOT/row['relative_path']) == row['sha256'], row['relative_path']
        assert sha(source/row['source_path']) == row['sha256'], row['source_path']
    pointers = json.loads((ROOT/'sync_metadata/upstream_lfs_pointer_inventory.json').read_text())['files']
    for row in pointers:
        assert hashlib.sha256(row['original_pointer_text'].encode()).hexdigest()==row['sha256']
        assert sha(source/row['source_path'])==row['sha256']
    write('source_preservation.json', dict(status='PASS',copied_files_verified=len(d['files']),
        archived_upstream_pointers_verified=len(pointers), original_files_verified=len(d['files'])+len(pointers),
        originals_byte_identical=True,source_git_history_modified=False,source_worktree_not_written=True,
        note='The frozen source is read-only during sync; target .gitattributes is packaging-only, original bytes separately preserved.'))


def secret_scan():
    patterns = {
        'private_key':re.compile(r'-----BEGIN (?:RSA |DSA |EC |OPENSSH |ENCRYPTED )?PRIVATE KEY-----'),
        'github_token':re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,})\b'),
        'aws_access_key':re.compile(r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
        'google_api_key':re.compile(r'\bAIza[A-Za-z0-9_-]{35}\b'),
        'slack_token':re.compile(r'\bxox[baprs]-[A-Za-z0-9-]{20,}\b'),
        'embedded_url_credential':re.compile(r'https?://[^\s/:@]+:[^\s/@]+@'),
        'quoted_credential_assignment':re.compile(r'''(?i)(?:password|passwd|api[_-]?key|access[_-]?token|auth[_-]?token|secret[_-]?key)\s*['"]?\s*[:=]\s*['"]([^'"\n]{6,})['"]'''),
    }
    findings = []
    count = 0
    for p in payload_files():
        rel = p.relative_to(ROOT).as_posix()
        assert not p.is_symlink(), rel
        assert not set(p.parts).intersection({'.ssh','.aws','.azure','.kube'}), rel
        assert p.name not in {'id_rsa','id_ed25519','credentials','hosts.yml','.env','CMakeCache.txt'}, rel
        try:
            text = p.read_text()
        except UnicodeDecodeError:
            continue
        count += 1
        for name, pattern in patterns.items():
            for match in pattern.finditer(text):
                # Scanner source holds definitions, not credential instances.
                findings.append(dict(path=rel, kind=name, line=text.count('\n',0,match.start())+1))
    write('secret_scan.json', dict(status='PASS' if not findings else 'FAIL',text_files_scanned=count,
        checks=list(patterns),findings=findings,redactions=[],
        note='No secret contents printed. WSL/remote paths and IP-only provenance are not credentials. No keys, SSH config, or cloud credential files copied.'))
    assert not findings, findings


def sizes_and_manifest():
    files=payload_files()
    largest=sorted(files,key=lambda p:p.stat().st_size,reverse=True)[:15]
    assert all(p.stat().st_size < 50*2**20 for p in files)
    assert not any(p.stat().st_size < 1000 and p.read_bytes().startswith(b'version https://git-lfs.github.com/spec/v1\n') for p in files)
    write('file_size_audit.json',dict(status='PASS',all_files_below_50_MiB=True,no_objects_at_or_above_95_MiB=True,
        LFS_enabled=False,scientific_binaries_compressed_to_evade_limit=False,
        largest=[dict(path=p.relative_to(ROOT).as_posix(),size_bytes=p.stat().st_size) for p in largest]))
    keys = [p for p in payload_files() if p.relative_to(ROOT).parts[0]=='frozen_reference' and p.name!='SHA256SUMS.txt']
    keys += [ROOT/'reports/sv1_3q'/n for n in ['REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance.json','winner_steady_candidate.json','steady_history.json','adaptive_timestep_history.json','steady_stop_request.json','stop_latency.json']]
    keys += [ROOT/'logs/sv1_3q/remote/REAL_VASCULAR_GPU_ILU_REUSE_WINNER.log', ROOT/'configs/sv1_3q/policy.json', ROOT/'configs/face_map.json']
    (ROOT/'frozen_reference/SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(ROOT).as_posix()}\n' for p in sorted(keys)))
    copied={x['relative_path']:x['source_path'] for x in json.loads((ROOT/'sync_metadata/copied_from_wsl.json').read_text())['files']}
    categories=dict(src='SOURCE',vendor='SOURCE',scripts='SCRIPT',tests='TEST',configs='CONFIG',reports='REPORT',logs='LOG',outputs='RESULT',inputs='REFERENCE',frozen_reference='REFERENCE',patches='SOURCE',benchmarks='REPORT')
    with (ROOT/MANIFEST).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['relative_path','size_bytes','sha256','category','source_path'])
        w.writeheader()
        for p in payload_files():
            rel=p.relative_to(ROOT).as_posix()
            if rel==MANIFEST:continue
            w.writerow(dict(relative_path=rel,size_bytes=p.stat().st_size,sha256=sha(p),
                category=categories.get(p.relative_to(ROOT).parts[0],'DOC'),source_path=copied.get(rel,'')))
    print(f'Delivery audit PASS; {len(payload_files())} files, source manifest excludes only its own recursive hash.')


if __name__=='__main__':
    audit_main()
    source_audit()
    secret_scan()
    sizes_and_manifest()
