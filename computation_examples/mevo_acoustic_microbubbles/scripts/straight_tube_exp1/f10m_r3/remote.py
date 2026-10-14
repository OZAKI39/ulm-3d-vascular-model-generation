"""Immutable, bounded supervisor execution; transfer only F10M_R3 inputs and source."""

import argparse
import io
import shlex
import subprocess
import tarfile
from .common import ROOT, BASE, identity, write

SERVER = "/workspace/mevo_acoustic_microbubbles/f10m_r3_dual_port"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("module", choices=["validate_gpu"])
    p.add_argument("--label", required=True)
    p.add_argument("--timeout", type=int, default=650)
    a = p.parse_args()
    sha = identity()["execution_commit"]
    paths = ["src", "configs/straight_tube_exp1", "scripts/straight_tube_exp1", "scripts/fast_pzt/python.sh"]
    dirty = subprocess.check_output(["git", "status", "--porcelain", "--", *paths], cwd=ROOT, text=True)
    if dirty:
        raise ValueError("Commit executable source before deployment: " + dirty)
    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    repo = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip()
    raw = subprocess.check_output(["git", "archive", sha, *[prefix + p for p in paths]], cwd=repo)
    buffer = io.BytesIO()
    with tarfile.open(fileobj=io.BytesIO(raw)) as src, tarfile.open(fileobj=buffer, mode="w:gz") as dst:
        for item in src:
            if item.name.startswith(prefix):
                item.name = item.name[len(prefix) :]
                if item.name:
                    dst.addfile(item, src.extractfile(item) if item.isfile() else None)
    archive = BASE / "execution" / (sha + ".tar.gz")
    archive.parent.mkdir(exist_ok=True)
    archive.write_bytes(buffer.getvalue())
    ex = SERVER + "/.execution/" + sha
    subprocess.run(
        ["ssh", "mevo-new", "mkdir -p " + shlex.quote(ex) + " && tar -xzf - -C " + shlex.quote(ex)],
        input=buffer.getvalue(),
        check=True,
    )
    # Only the small response table is needed; no Y field, mesh or CFD dependency.
    subprocess.run(
        ["ssh", "mevo-new", "mkdir -p " + SERVER + "/results/straight_tube_exp1/f10m_r3/inputs"], check=True
    )
    subprocess.run(
        [
            "rsync",
            "-a",
            str(BASE / "inputs/response.h5"),
            "mevo-new:" + SERVER + "/results/straight_tube_exp1/f10m_r3/inputs/",
        ],
        check=True,
    )
    for rel in [
        "CONFIG_SNAPSHOT.json",
        "INPUT_IDENTITY.json",
        "initial_state.npz",
        "wall_reference.npz",
        "analysis/ACOUSTIC_ADMISSION.json",
    ]:
        subprocess.run(
            [
                "rsync",
                "-a",
                "--relative",
                str(BASE) + "/./" + rel,
                "mevo-new:" + SERVER + "/results/straight_tube_exp1/f10m_r3/",
            ],
            check=True,
        )
    name = "mevo_f10m_r3_" + a.label
    command = shlex.join(
        [
            "/usr/bin/timeout",
            str(a.timeout),
            "/bin/bash",
            "scripts/fast_pzt/python.sh",
            "-m",
            "scripts.straight_tube_exp1.f10m_r3." + a.module,
        ]
    )
    data = dict(ex=ex, server=SERVER, sha=sha, name=name, label=a.label, command=command)
    code = """from pathlib import Path
import json,subprocess
d=DATA
ex=Path(d['ex']);root=Path(d['server'])
assert not subprocess.check_output(['nvidia-smi','--query-compute-apps=pid','--format=csv,noheader'],text=True).strip(),'GPU occupied; do not interrupt other work'
for name,target in [('results',root/'results'),('.venv',root.parent/'.venv')]:
 p=ex/name
 if not p.exists(): p.symlink_to(target,target_is_directory=True)
(ex/'SOURCE_F10M_R3_COMMIT.txt').write_text(d['sha']+'\\n')
out=root/'results/straight_tube_exp1/f10m_r3/logs';out.mkdir(exist_ok=True,parents=True)
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf')
assert not conf.exists(),'Use a new run label; preserve failed attempts'
conf.write_text('[program:'+d['name']+']\\ncommand='+d['command']+'\\ndirectory='+str(ex)+'\\nautostart=false\\nautorestart=false\\nstartsecs=0\\nstopasgroup=true\\nkillasgroup=true\\nstdout_logfile='+str(out/(d['label']+'.log'))+'\\nstderr_logfile='+str(out/(d['label']+'.stderr.log'))+'\\nstdout_logfile_maxbytes=0\\nstderr_logfile_maxbytes=0\\nenvironment=F10M_R3_ATTEMPT="'+d['label']+'",PYTHONUNBUFFERED="1",WARP_CACHE_PATH="/workspace/mevo_acoustic_microbubbles/.cache/warp",XDG_CACHE_HOME="/workspace/mevo_acoustic_microbubbles/.cache"\\n')
for args in [['reread'],['update',d['name']],['start',d['name']]]:subprocess.run(['supervisorctl',*args],check=True)
(out/(d['label']+'_EXECUTION.json')).write_text(json.dumps(d,indent=2))
print(json.dumps(d))
""".replace("DATA", repr(data))
    subprocess.run(["ssh", "mevo-new", "python3 -"], input=code, text=True, check=True)
    write(BASE / "execution" / (a.label + ".json"), data)


if __name__ == "__main__":
    main()
