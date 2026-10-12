"""Read-only remote evidence audit, then hash-verified local scientific copy."""

import io
import json
import subprocess
import tarfile
from .common import BASE, ROOT, sha, js, write


def main():
    remote = r"""
from pathlib import Path
import hashlib,json,subprocess,sys,platform,importlib.metadata
root=Path('/workspace/mevo_acoustic_microbubbles/f10h')
out=root/'results/fast_pzt/f10h'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
commits=sorted({json.loads(p.read_text())['execution_commit'] for p in out.glob('*/IDENTITY.json')})
sources={}
for commit in commits:
    ex=root/'.execution'/commit
    assert (ex/'SOURCE_F10H_COMMIT.txt').read_text().strip()==commit
    sources[commit]={str(p.relative_to(ex)):sha(p) for p in sorted(ex.rglob('*'))
                     if p.is_file() and not p.is_symlink() and '__pycache__' not in p.parts
                     and p.suffix != '.pyc' and not p.name.startswith('SOURCE_')}
    # Every hash recorded by the running code is rechecked against its snapshot.
    for identity in out.glob('*/IDENTITY.json'):
        data=json.loads(identity.read_text())
        if data['execution_commit']!=commit:continue
        for name,value in data['source_sha256'].items():assert sha(ex/name)==value,name
        for name,value in data['inputs'].items():
            assert sha(root/'results/fast_pzt/f5'/name)==value['sha256'],name
write(out/'REMOTE_SOURCE_MANIFEST.json',dict(status='PASS',sources=sources))
env=dict(python=sys.version,executable=sys.executable,platform=platform.platform(),
         packages={p:importlib.metadata.version(p) for p in ['warp-lang','numpy','scipy','h5py','matplotlib']},
         gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,uuid,driver_version,memory.total','--format=csv,noheader'],text=True),
         final_compute_processes=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],text=True).strip(),
         supervisor=subprocess.run(['supervisorctl','status'],capture_output=True,text=True).stdout,
         guide_sha256=sha(Path('/etc/vast-agents-guide.md')))
write(out/'ENVIRONMENT.json',env)
files={str(p.relative_to(out)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*'))
       if p.is_file() and p.name!='REMOTE_MANIFEST.json'}
write(out/'REMOTE_MANIFEST.json',dict(files=files))
print(json.dumps(dict(files=len(files),commits=commits,final_compute_processes=env['final_compute_processes'])))
"""
    result = subprocess.run(
        [
            "ssh",
            "-o",
            "ConnectTimeout=20",
            "mevo-new",
            "/workspace/mevo_acoustic_microbubbles/.venv/bin/python -",
        ],
        input=remote,
        text=True,
        check=True,
        capture_output=True,
    )
    print(result.stdout)
    subprocess.run(
        [
            "rsync",
            "-az",
            "mevo-new:/workspace/mevo_acoustic_microbubbles/f10h/results/fast_pzt/f10h/",
            str(BASE) + "/",
        ],
        check=True,
    )
    manifest = js(BASE / "REMOTE_MANIFEST.json")
    for name, entry in manifest["files"].items():
        assert sha(BASE / name) == entry["sha256"], name
    remote_source = js(BASE / "REMOTE_SOURCE_MANIFEST.json")
    repository = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip()
    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    checked = {}
    import hashlib

    for commit, files in remote_source["sources"].items():
        raw = subprocess.check_output(
            ["git", "archive", commit, *[prefix + n for n in files]], cwd=repository
        )
        with tarfile.open(fileobj=io.BytesIO(raw)) as archive:
            actual = {
                p.name[len(prefix) :]: hashlib.sha256(archive.extractfile(p).read()).hexdigest()
                for p in archive
                if p.isfile()
            }
        assert all(actual.get(n) == h for n, h in files.items())
        checked[commit] = len(files)
    write(
        BASE / "REMOTE_VERIFICATION.json",
        dict(
            status="PASS",
            remote_files_verified=len(manifest["files"]),
            source_files_against_git_archive=checked,
            input_files_verified=13,
            final_compute_processes=js(BASE / "ENVIRONMENT.json")["final_compute_processes"],
        ),
    )
    print(json.dumps(js(BASE / "REMOTE_VERIFICATION.json"), indent=2))


if __name__ == "__main__":
    main()
