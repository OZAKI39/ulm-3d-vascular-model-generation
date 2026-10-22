"""Commit-only deployment to an isolated F10J root and unique supervisor jobs."""

import argparse
import io
from pathlib import Path
import shlex
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[2]
SERVER = "/workspace/mevo_acoustic_microbubbles/f10j"
PREFIX = "computation_examples/mevo_acoustic_microbubbles/"
SSH = ["ssh", "-o", "ConnectTimeout=20", "-o", "ServerAliveInterval=15", "mevo-new"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("module")
    ap.add_argument("--label", required=True)
    ap.add_argument("--timeout", type=int, default=240)
    ap.add_argument(
        "--deploy-only", action="store_true", help="Install immutable committed source without starting a job"
    )
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
            "scripts/fast_pzt_f6",
            "scripts/fast_pzt_f7",
            "scripts/fast_pzt_f8",
            "scripts/fast_pzt_f10a",
            "scripts/fast_pzt_f10b",
            "scripts/fast_pzt_f10c",
            "scripts/fast_pzt_f10d",
            "scripts/fast_pzt_f10e",
            "scripts/fast_pzt_f10f",
            "scripts/fast_pzt_f10g",
            "scripts/fast_pzt_f10h",
            "scripts/fast_pzt_f10i",
            "scripts/fast_pzt_f10j",
            "configs/fast_pzt/f10j",
        ],
        cwd=ROOT,
        text=True,
    )
    if dirty:
        raise ValueError("Commit executable F10J source before deployment: " + dirty)
    repo = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip()
    paths = [
        "src",
        "scripts/fast_pzt",
        "scripts/fast_pzt_f2",
        "scripts/fast_pzt_f3",
        "scripts/fast_pzt_f4",
        "scripts/fast_pzt_f5",
        "scripts/fast_pzt_f6",
        "scripts/fast_pzt_f7",
        "scripts/fast_pzt_f8",
        "scripts/fast_pzt_f9",
        "scripts/fast_pzt_f10a",
        "scripts/fast_pzt_f10b",
        "scripts/fast_pzt_f10c",
        "scripts/fast_pzt_f10d",
        "scripts/fast_pzt_f10e",
        "scripts/fast_pzt_f10f",
        "scripts/fast_pzt_f10g",
        "scripts/fast_pzt_f10h",
        "scripts/fast_pzt_f10i",
        "scripts/fast_pzt_f10j",
        "validation/history_force_sensitivity",
        "configs",
        "docs/FAST_PZT_F10A_DEVELOPMENT.md",
        "docs/FAST_PZT_F10C_DEVELOPMENT.md",
        "docs/FAST_PZT_F10C_PERFORMANCE.md",
        "docs/FAST_PZT_F10D_DEVELOPMENT.md",
        "docs/FAST_PZT_F10D_PERFORMANCE.md",
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
    marker = execution + "/SOURCE_F10J_COMMIT.txt"
    exists = subprocess.run(SSH + ["test -f " + shlex.quote(marker)]).returncode == 0
    if not exists:
        subprocess.run(
            SSH + ["mkdir -p " + shlex.quote(execution) + " && tar -xzf - -C " + shlex.quote(execution)],
            input=data.getvalue(),
            check=True,
        )
    name = "mevo_fast_pzt_f10j_" + a.label
    command = shlex.join(
        [
            "/usr/bin/timeout",
            str(a.timeout),
            "/bin/bash",
            "scripts/fast_pzt/python.sh",
            "-m",
            "scripts.fast_pzt_f10j.timed",
            "--module",
            a.module,
            "--label",
            a.label,
            *extra,
        ]
    )
    d = dict(
        root=SERVER,
        execution=execution,
        sha=sha,
        name=name,
        label=a.label,
        command=command,
        deploy_only=a.deploy_only,
    )
    code = """from pathlib import Path
import subprocess,json
d=DATA
root=Path(d['root']);ex=Path(d['execution'])
for target in ['results','.venv']:
 p=ex/target
 if not p.exists(): p.symlink_to(root/target if target=='results' else root.parent/'.venv',target_is_directory=True)
marker=ex/'SOURCE_F10J_COMMIT.txt'
if marker.exists(): assert marker.read_text().strip()==d['sha']
else: marker.write_text(d['sha']+'\\n')
(ex/'SOURCE_F10C_COMMIT.txt').write_text(d['sha']+'\\n')
if d['deploy_only']:
 print(json.dumps(d))
 raise SystemExit(0)
states=subprocess.run(['supervisorctl','status'],capture_output=True,text=True).stdout
assert not any(line.startswith('mevo_fast_') and 'RUNNING' in line for line in states.splitlines()), 'Another MEVO GPU science task is running'
gpu_processes=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid','--format=csv,noheader'],text=True).strip()
assert not gpu_processes, 'GPU is occupied by another process; wait for it to finish, do not interrupt unrelated jobs: '+gpu_processes
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf')
assert not conf.exists(),'New run label required; failed attempts are retained'
out=root/'results/fast_pzt/f10j'
out.mkdir(parents=True,exist_ok=True)
for stage in ['f1','f2','f3','f4','f5','f6','f7','f8','f9','f10a','f10b','f10c','f10d','f10e','f10f','f10g','f10h','f10i']:
 p=root/'results/fast_pzt'/stage
 if not p.exists(): p.symlink_to(root.parent/(stage+'/results/fast_pzt/'+stage if stage in ['f6','f7','f8','f9','f10a','f10b','f10c','f10d','f10e','f10f','f10g','f10h','f10i'] else 'project/results/fast_pzt/'+stage),target_is_directory=True)
conf.write_text('[program:'+d['name']+']\\ncommand='+d['command']+'\\ndirectory='+str(ex)+'\\nautostart=false\\nautorestart=false\\nstartsecs=0\\nstopasgroup=true\\nkillasgroup=true\\nstdout_logfile='+str(out/(d['label']+'.log'))+'\\nstderr_logfile='+str(out/(d['label']+'.stderr.log'))+'\\nstdout_logfile_maxbytes=0\\nstderr_logfile_maxbytes=0\\nenvironment=PYTHONUNBUFFERED="1",WARP_CACHE_PATH="/workspace/mevo_acoustic_microbubbles/.cache/warp",XDG_CACHE_HOME="/workspace/mevo_acoustic_microbubbles/.cache"\\n')
for args in [['reread'],['update',d['name']],['start',d['name']]]: subprocess.run(['supervisorctl',*args],check=True)
print(json.dumps(d))
""".replace("DATA", repr(d))
    subprocess.run(SSH + ["python3 -"], input=code, text=True, check=True)


if __name__ == "__main__":
    main()
