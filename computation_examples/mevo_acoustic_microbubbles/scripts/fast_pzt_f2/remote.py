"""Deploy committed source only; run bounded F1 jobs under server supervisor."""

import argparse
from pathlib import Path
import shlex
import subprocess
import tarfile
import io

ROOT = Path(__file__).resolve().parents[2]
SERVER = "/workspace/mevo_fast_pzt_f2"
PREFIX = "computation_examples/mevo_acoustic_microbubbles/"
SSH = [
    "ssh",
    "-o",
    "ConnectTimeout=20",
    "-o",
    "ServerAliveInterval=15",
    "-o",
    "ServerAliveCountMax=3",
    "vast4090",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("module")
    ap.add_argument("--label", required=True)
    ap.add_argument("--timeout", type=int, default=240)
    a, extra = ap.parse_known_args()
    if any(c not in "abcdefghijklmnopqrstuvwxyz0123456789_" for c in a.label):
        raise ValueError("Safe unique label required")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = subprocess.check_output(
        ["git", "diff", "--name-only", "HEAD", "--", "src/mevo_mb/fast_sim", "scripts/fast_pzt_f2"],
        cwd=ROOT,
        text=True,
    )
    if dirty:
        raise ValueError("Commit F1 executable sources before deployment: " + dirty)
    repo = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip()
    raw = subprocess.check_output(
        [
            "git",
            "archive",
            head,
            *[
                PREFIX + p
                for p in [
                    "src",
                    "scripts/fast_pzt",
                    "scripts/fast_pzt_f2",
                    "configs",
                    "docs/FAST_PZT_F2_PLAN.md",
                ]
            ],
        ],
        cwd=repo,
    )
    packed = io.BytesIO()
    with tarfile.open(fileobj=io.BytesIO(raw)) as source, tarfile.open(fileobj=packed, mode="w:gz") as target:
        for member in source:
            if member.name.startswith(PREFIX):
                member.name = member.name[len(PREFIX) :]
                if member.name:
                    target.addfile(member, source.extractfile(member) if member.isfile() else None)
    execution = SERVER + "/.execution/" + head
    command = "mkdir -p " + shlex.quote(execution) + " && tar -xzf - -C " + shlex.quote(execution)
    # Deployment targets are immutable after first extraction.
    exists = (
        subprocess.run(SSH + ["test -f " + shlex.quote(execution + "/SOURCE_F2_COMMIT.txt")]).returncode == 0
    )
    if not exists:
        subprocess.run(SSH + [command], input=packed.getvalue(), check=True)
    name = "mevo_fast_pzt_f2_" + a.label
    extra = extra[1:] if extra[:1] == ["--"] else extra
    args = [
        "/usr/bin/timeout",
        str(a.timeout),
        "/bin/bash",
        "scripts/fast_pzt/python.sh",
        "-m",
        "scripts.fast_pzt_f2." + a.module,
        "--label",
        a.label,
        *extra,
    ]
    script = """
from pathlib import Path
import subprocess,json
d=DATA
root=Path(d['root']); execution=Path(d['execution'])
for target in ['results','.venv']:
 p=execution/target
 if not p.exists(): p.symlink_to(root/target,target_is_directory=True)
p=execution/'SOURCE_F2_COMMIT.txt'
if p.exists(): assert p.read_text().strip()==d['head']
else: p.write_text(d['head']+'\\n')
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf'); assert not conf.exists(),'Use a new run label'
out=root/'results/fast_pzt/f2'
conf.write_text('[program:'+d['name']+']\\ncommand='+d['command']+'\\ndirectory='+str(execution)+'\\nautostart=false\\nautorestart=false\\nstartsecs=0\\nstopasgroup=true\\nkillasgroup=true\\nstdout_logfile='+str(out/(d['label']+'.log'))+'\\nstderr_logfile='+str(out/(d['label']+'.stderr.log'))+'\\nstdout_logfile_maxbytes=0\\nstderr_logfile_maxbytes=0\\nenvironment=PYTHONUNBUFFERED="1"\\n')
for args in [['reread'],['update',d['name']],['start',d['name']]]: subprocess.run(['supervisorctl',*args],check=True)
print(json.dumps(d))
"""
    data = dict(
        root=SERVER, execution=execution, head=head, name=name, label=a.label, command=shlex.join(args)
    )
    script = script.replace("DATA", repr(data))
    subprocess.run(SSH + ["python3 -"], input=script, text=True, check=True)


if __name__ == "__main__":
    main()
