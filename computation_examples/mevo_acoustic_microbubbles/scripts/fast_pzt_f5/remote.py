"""Commit-only deployment to an isolated F5 root and unique supervisor jobs."""

import argparse
import io
from pathlib import Path
import shlex
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[2]
SERVER = "/workspace/mevo_fast_pzt_f5"
PREFIX = "computation_examples/mevo_acoustic_microbubbles/"
SSH = ["ssh", "-o", "ConnectTimeout=20", "-o", "ServerAliveInterval=15", "vast4090"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("module")
    ap.add_argument("--label", required=True)
    ap.add_argument("--timeout", type=int, default=240)
    a, extra = ap.parse_known_args()
    if any(c not in "abcdefghijklmnopqrstuvwxyz0123456789_" for c in a.label + a.module):
        raise ValueError("Safe unique module/label required")
    sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = subprocess.check_output(
        [
            "git",
            "status",
            "--porcelain",
            "--",
            "src/mevo_mb/fast_sim",
            "scripts/fast_pzt_f5",
            "configs/fast_pzt/f5",
        ],
        cwd=ROOT,
        text=True,
    )
    if dirty:
        raise ValueError("Commit executable F5 source before deployment: " + dirty)
    repo = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip()
    paths = [
        "src",
        "scripts/fast_pzt",
        "scripts/fast_pzt_f2",
        "scripts/fast_pzt_f3",
        "scripts/fast_pzt_f4",
        "scripts/fast_pzt_f5",
        "configs",
        "docs/FAST_PZT_F5_PLAN.md",
    ]
    raw = subprocess.check_output(["git", "archive", sha, *[PREFIX + p for p in paths]], cwd=repo)
    data = io.BytesIO()
    with tarfile.open(fileobj=io.BytesIO(raw)) as source, tarfile.open(fileobj=data, mode="w:gz") as target:
        for item in source:
            if item.name.startswith(PREFIX):
                item.name = item.name[len(PREFIX) :]
                if item.name:
                    target.addfile(item, source.extractfile(item) if item.isfile() else None)
    execution = SERVER + "/.execution/" + sha
    marker = execution + "/SOURCE_F5_COMMIT.txt"
    exists = subprocess.run(SSH + ["test -f " + shlex.quote(marker)]).returncode == 0
    if not exists:
        subprocess.run(
            SSH + ["mkdir -p " + shlex.quote(execution) + " && tar -xzf - -C " + shlex.quote(execution)],
            input=data.getvalue(),
            check=True,
        )
    name = "mevo_fast_pzt_f5_" + a.label
    command = shlex.join(
        [
            "/usr/bin/timeout",
            str(a.timeout),
            "/bin/bash",
            "scripts/fast_pzt/python.sh",
            "-m",
            "scripts.fast_pzt_f5." + a.module,
            "--label",
            a.label,
            *extra,
        ]
    )
    d = dict(root=SERVER, execution=execution, sha=sha, name=name, label=a.label, command=command)
    code = """from pathlib import Path
import subprocess,json
d=DATA
root=Path(d['root']);ex=Path(d['execution'])
for target in ['results','.venv']:
 p=ex/target
 if not p.exists(): p.symlink_to(root/target,target_is_directory=True)
marker=ex/'SOURCE_F5_COMMIT.txt'
if marker.exists(): assert marker.read_text().strip()==d['sha']
else: marker.write_text(d['sha']+'\\n')
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf')
assert not conf.exists(),'New run label required; failed attempts are retained'
out=root/'results/fast_pzt/f5'
conf.write_text('[program:'+d['name']+']\\ncommand='+d['command']+'\\ndirectory='+str(ex)+'\\nautostart=false\\nautorestart=false\\nstartsecs=0\\nstopasgroup=true\\nkillasgroup=true\\nstdout_logfile='+str(out/(d['label']+'.log'))+'\\nstderr_logfile='+str(out/(d['label']+'.stderr.log'))+'\\nstdout_logfile_maxbytes=0\\nstderr_logfile_maxbytes=0\\nenvironment=PYTHONUNBUFFERED="1"\\n')
for args in [['reread'],['update',d['name']],['start',d['name']]]: subprocess.run(['supervisorctl',*args],check=True)
print(json.dumps(d))
""".replace("DATA", repr(d))
    subprocess.run(SSH + ["python3 -"], input=code, text=True, check=True)


if __name__ == "__main__":
    main()
