"""Mirror this new stage only; verify assets and unchanged remote donors, without running physics."""
from pathlib import Path, PurePosixPath
import argparse, datetime, hashlib, json, shlex, subprocess

S = Path(__file__).resolve().parents[1]
arg = argparse.ArgumentParser()
arg.add_argument('--ssh-alias', required=True)
args = arg.parse_args()
paths = json.loads((S / 'provenance/STAGE_PATHS.json').read_text())
remote = paths['remote_stage']
assert PurePosixPath(remote).name == 'microbubble_vascular_workflow_v0_' + S.name
assert str(S) == paths['local_stage']
excluded = {'provenance/SYNC_MANIFEST.json', 'validation/REMOTE_LOCAL_SHA256_CHECK.json', 'validation/SYNC_TRANSPORT.log'}

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

files = {}
for p in sorted(S.rglob('*')):
    assert not p.is_symlink(), 'Unexpected symlink: ' + str(p)
    if p.is_file() and str(p.relative_to(S)) not in excluded:
        assert p.suffix != '.pyc' and '__pycache__' not in p.parts
        files[str(p.relative_to(S))] = {'sha256': digest(p), 'bytes': p.stat().st_size}
manifest = {'scope': 'All new-stage regular files except named transport evidence files', 'excluded': sorted(excluded), 'files': files}
m = S / 'provenance/SYNC_MANIFEST.json'
m.write_text(json.dumps(manifest, indent=2) + '\n')
rsync = subprocess.run(['rsync', '-az', '--exclude=__pycache__', '--exclude=*.pyc', str(S) + '/', args.ssh_alias + ':' + remote + '/'], text=True, capture_output=True, check=True)
log = rsync.stdout + rsync.stderr
probe = json.loads((S / 'provenance/REMOTE_INPUTS_TO_PRESERVE.json').read_text())
identity = json.loads((S / 'provenance/BASELINE_CPU_IDENTITY_TO_RECHECK.json').read_text())
payload = {'remote': remote, 'manifest_sha256': digest(m), 'original_inputs': probe['inputs'], 'original_builds': identity['binaries'], 'hostname': probe['hostname']}
code = '''from pathlib import Path
import json,hashlib,platform,datetime
payload = PAYLOAD
S=Path(payload['remote'])
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
m=S/'provenance/SYNC_MANIFEST.json'
assert digest(m)==payload['manifest_sha256']
manifest=json.loads(m.read_text()); errors=[]
for rel,d in manifest['files'].items():
 p=S/rel
 if not p.is_file() or digest(p)!=d['sha256'] or p.stat().st_size!=d['bytes']: errors.append(rel)
actual={str(p.relative_to(S)) for p in S.rglob('*') if p.is_file() and str(p.relative_to(S)) not in manifest['excluded']}
missing=set(manifest['files'])-actual;extra=actual-set(manifest['files'])
changed=[]
for path,d in payload['original_inputs'].items():
 if digest(Path(path))!=d['sha256']:changed.append(path)
for name,d in payload['original_builds'].items():
 if digest(Path(d['path']))!=d['sha256']:changed.append(d['path'])
host_ok=platform.node()==payload['hostname']
result={'REMOTE_LOCAL_SHA256_STATUS':'PASS' if not(errors or missing or extra or changed) and host_ok else 'FAIL','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'hostname':platform.node(),'hostname_unchanged':host_ok,'remote_stage':str(S),'manifest_sha256':digest(m),'checked_stage_files':len(manifest['files']),'checked_stage_bytes':sum(x['bytes'] for x in manifest['files'].values()),'mismatched_files':errors,'missing_files':sorted(missing),'unexpected_files':sorted(extra),'original_remote_inputs_rechecked':len(payload['original_inputs']),'original_CPU_build_files_rechecked':len(payload['original_builds']),'changed_original_remote_files':changed,'excluded_transport_evidence':manifest['excluded'],'manifest_verified_separately':True,'receipt_transport_note':'This receipt is subsequently copied to the remote new stage and its SHA256 compared; terminal/transport log records that final comparison.'}
print(json.dumps(result,indent=2))
'''.replace('PAYLOAD', repr(payload))
r = subprocess.run(['ssh', args.ssh_alias, 'python3', '-'], input=code, text=True, capture_output=True, check=True)
log += r.stderr
result = json.loads(r.stdout)
assert result['REMOTE_LOCAL_SHA256_STATUS'] == 'PASS', result
# Verify local assets did not change during transfer.
assert all(digest(S / rel) == d['sha256'] for rel, d in files.items())
result['local_manifest_rechecked_after_transfer'] = True
receipt = S / 'validation/REMOTE_LOCAL_SHA256_CHECK.json'
receipt.write_text(json.dumps(result, indent=2) + '\n')
c = subprocess.run(['scp', str(receipt), args.ssh_alias + ':' + remote + '/validation/'], text=True, capture_output=True, check=True)
log += c.stderr
r = subprocess.run(['ssh', args.ssh_alias, 'sha256sum ' + shlex.quote(remote + '/validation/REMOTE_LOCAL_SHA256_CHECK.json')], text=True, capture_output=True, check=True)
assert r.stdout.split()[0] == digest(receipt)
log += r.stderr + r.stdout + 'RECEIPT_REMOTE_LOCAL_SHA256_STATUS: PASS\n'
(S / 'validation/SYNC_TRANSPORT.log').write_text(log)
c = subprocess.run(['scp', str(S / 'validation/SYNC_TRANSPORT.log'), args.ssh_alias + ':' + remote + '/validation/'], text=True, capture_output=True, check=True)
r = subprocess.run(['ssh', args.ssh_alias, 'sha256sum ' + shlex.quote(remote + '/validation/SYNC_TRANSPORT.log')], text=True, capture_output=True, check=True)
assert r.stdout.split()[0] == digest(S / 'validation/SYNC_TRANSPORT.log')
print(json.dumps(result, indent=2))
print('RECEIPT_REMOTE_LOCAL_SHA256_STATUS: PASS')
print('TRANSPORT_LOG_REMOTE_LOCAL_SHA256_STATUS: PASS')
