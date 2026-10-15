"""Deploy committed source and hashed inputs to an immutable supervised execution."""

import argparse
import io
import shlex
import subprocess
import tarfile
import gzip
from .common import ROOT, BASE, OLD, CONFIG, identity, read, write, sha

SERVER = "/workspace/mevo_acoustic_microbubbles/f10s_control_validation"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--label", required=True)
    p.add_argument("--mode", choices=["A","B"], default="A")
    p.add_argument("--timeout", type=int, default=340)
    p.add_argument("--route", choices=["A", "B", "reserve"], default="A")
    a = p.parse_args()
    commit = identity()["execution_commit"]
    paths = ["src", "configs/straight_tube_exp1", "scripts/straight_tube_exp1", "scripts/fast_pzt/python.sh"]
    assert not subprocess.check_output(
        ["git", "status", "--porcelain", "--", *paths], cwd=ROOT, text=True
    ).strip()
    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    repo = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip()
    raw = subprocess.check_output(["git", "archive", commit, *[prefix + k for k in paths]], cwd=repo)
    buffer = io.BytesIO()
    with tarfile.open(fileobj=io.BytesIO(raw)) as src, tarfile.open(fileobj=buffer, mode="w:gz") as dst:
        for item in src:
            if item.name.startswith(prefix):
                item.name = item.name[len(prefix) :]
                if item.name:
                    dst.addfile(item, src.extractfile(item) if item.isfile() else None)
    archive = BASE / "execution" / f"{commit}.tar.gz"
    archive.parent.mkdir(exist_ok=True)
    if archive.exists():
        assert gzip.decompress(archive.read_bytes()) == gzip.decompress(buffer.getvalue())
        buffer = io.BytesIO(archive.read_bytes())
    else:
        archive.write_bytes(buffer.getvalue())
    ex = SERVER + "/.execution/" + commit
    subprocess.run(
        ["ssh", "mevo-new", "mkdir -p " + shlex.quote(ex) + " && tar -xzf - -C " + shlex.quote(ex)],
        input=buffer.getvalue(),
        check=True,
    )
    init = """from pathlib import Path
root=Path(SERVER); old=Path('/workspace/mevo_acoustic_microbubbles/f10r_scheduled_control')
data=root/'results/straight_tube_exp1';data.mkdir(parents=True,exist_ok=True)
for name in ['f10n','f10q','f10r']:
 target=old/'results/straight_tube_exp1'/name
 assert target.exists(),str(target)
 link=data/name
 if not link.exists(): link.symlink_to(target,target_is_directory=True)
(data/'f10s').mkdir(exist_ok=True)
""".replace('SERVER',repr(SERVER))
    subprocess.run(['ssh','mevo-new','python3 -'],input=init,text=True,check=True)
    for name in ['PREREGISTRATION.json','INITIAL_STATES.npz','RUN_SCHEDULES.json']:
        subprocess.run(['rsync','-a',str(BASE/name),'mevo-new:'+SERVER+'/results/straight_tube_exp1/f10s/'],check=True)
    label = "mevo_f10s_" + a.label
    command = shlex.join(
        [
            "/usr/bin/timeout",
            "--signal=TERM",
            "--kill-after=10",
            str(a.timeout),
            "/bin/bash",
            "scripts/fast_pzt/python.sh",
            "-m",
            "scripts.straight_tube_exp1.f10s.run",
        ]
    )
    data = dict(
        ex=ex,
        server=SERVER,
        sha=commit,
        name=label,
        label=a.label,
        mode=a.mode,
        route=a.mode,
        command=command,
        archive_sha256=sha(archive),
        protocol_sha256=sha(CONFIG),
    )
    code = """from pathlib import Path
import subprocess,json
d=DATA
ex=Path(d['ex']); root=Path(d['server'])
assert not subprocess.check_output(['nvidia-smi','--query-compute-apps=pid','--format=csv,noheader'],text=True).strip(),'GPU occupied; no other job will be stopped'
for name,target in [('results',root/'results'),('.venv',root.parent/'.venv')]:
 p=ex/name
 if not p.exists():p.symlink_to(target,target_is_directory=True)
(ex/'SOURCE_F10S_COMMIT.txt').write_text(d['sha']+'\\n')
out=root/'results/straight_tube_exp1/f10s/logs';out.mkdir(parents=True,exist_ok=True)
conf=Path('/etc/supervisor/conf.d')/(d['name']+'.conf')
assert not conf.exists(),'New immutable label required'
conf.write_text('[program:'+d['name']+']\\ncommand='+d['command']+'\\ndirectory='+str(ex)+'\\nautostart=false\\nautorestart=false\\nstartsecs=0\\nstopasgroup=true\\nkillasgroup=true\\nstdout_logfile='+str(out/(d['label']+'.log'))+'\\nstderr_logfile='+str(out/(d['label']+'.stderr.log'))+'\\nstdout_logfile_maxbytes=0\\nstderr_logfile_maxbytes=0\\nenvironment=F10S_ROUTE="'+d['route']+'",F10S_ATTEMPT="'+d['label']+'",F10S_MODE="'+d['mode']+'",PYTHONUNBUFFERED="1",WARP_CACHE_PATH="/workspace/mevo_acoustic_microbubbles/.cache/warp",XDG_CACHE_HOME="/workspace/mevo_acoustic_microbubbles/.cache"\\n')
(out/(d['label']+'_supervisor.conf')).write_text(conf.read_text())
for args in [['reread'],['update',d['name']],['start',d['name']]]:subprocess.run(['supervisorctl',*args],check=True)
(out/(d['label']+'_EXECUTION.json')).write_text(json.dumps(d,indent=2))
print(json.dumps(d))
""".replace("DATA", repr(data))
    subprocess.run(["ssh", "mevo-new", "python3 -"], input=code, text=True, check=True)
    write(BASE / "execution" / f"{a.label}.json", data)


if __name__ == "__main__":
    main()
