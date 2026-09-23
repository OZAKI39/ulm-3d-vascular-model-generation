#!/usr/bin/env python3
"""WSL-owned remote execution. No remote editing and no rsync deletion."""
import argparse
import json
import os
import platform
import resource
import shlex
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from fem3d.audit import command, sha256, timestamp, write_json

OWNER = "formal-fem3d-wsl-source-of-truth-stage00"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=("probe", "sync", "run", "fetch"))
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--stage", choices=("0", "1", "1.5", "1.6", "1.7", "1.8", "2", "3"), default="0")
    p.add_argument("--label", default="remote_command")
    p.add_argument("--ranks", type=int, default=1)
    p.add_argument("--python", default="python3", help="Remote interpreter for provenance wrapper")
    p.add_argument("argv", nargs="*")
    # Options after the action remain supported; command is explicitly after --.
    raw = sys.argv[1:]
    sep = raw.index("--") if "--" in raw else len(raw)
    a = p.parse_args(raw[:sep])
    stage = {"1.5": "stage01_5", "1.6": "stage01_6", "1.7": "stage01_7", "1.8": "stage01_8"}[a.stage] if a.stage in ("1.5", "1.6", "1.7", "1.8") else f"stage{int(a.stage):02d}"
    user_command = raw[sep+1:]
    cfg = json.loads((ROOT/"configs/remote.json").read_text())
    if ROOT != Path(cfg["source_of_truth"]).resolve():
        raise RuntimeError("Remote orchestration must run from the configured WSL source-of-truth directory")
    ssh = ["ssh", "-p", str(cfg["port"]), "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
           "-o", "ServerAliveInterval=30", "-o", "ServerAliveCountMax=6",
           "-o", "IdentitiesOnly=yes", "-i", cfg["identity_file"]]
    state_path = ROOT/"remote/connection.local.json"
    log = ROOT/"logs"/stage/(time.strftime("%Y%m%dT%H%M%S", time.gmtime())+f"_{time.time_ns()%1000000000}_remote_{a.action}")
    log.mkdir(parents=True)
    evidence = {"command": sys.argv, "timestamp": timestamp(), "host": cfg["host"],
                "hostname":platform.node(), "python_version":sys.version,
                "git":command(["git","-C",str(ROOT),"rev-parse","HEAD"]),
                "configuration_snapshot": cfg, "dry_run": a.dry_run, "commands": [],
                "mpi_ranks": a.ranks, "script_sha256": sha256(__file__)}
    for name in ("local_audit_environment", "remote_environment", "remote_fem_environment"):
        snapshot=ROOT/f"reports/stage00/{name}.json"
        if snapshot.exists():
            evidence[name+"_snapshot"] = {"path":str(snapshot),"sha256":sha256(snapshot),
                                          "data":json.loads(snapshot.read_text())}
    start = time.perf_counter()

    def execute(cmd, stdin=None):
        print("$", shlex.join(cmd), flush=True)
        t = time.perf_counter()
        r = subprocess.run(cmd, input=stdin, text=True, capture_output=True)
        evidence["commands"].append({"command": cmd, "returncode": r.returncode,
                                     "wall_time_s": time.perf_counter()-t})
        with (log/"stdout.log").open("a") as f:
            f.write(r.stdout)
        with (log/"stderr.log").open("a") as f:
            f.write(r.stderr)
        if r.returncode:
            print(r.stdout[-4000:]); print(r.stderr[-4000:], file=sys.stderr)
            raise subprocess.CalledProcessError(r.returncode, cmd)
        return r.stdout

    def ssh_python(code):
        return execute(ssh+[cfg["host"], "python3 -B -"], code)

    try:
        if a.action == "probe":
            audit = (ROOT/"src/fem3d/audit.py").read_text()
            code = audit + "\nimport shutil, glob\n" + f"config={cfg!r}\n"
            code += """
data = environment()
candidates = []
for base in (config['preferred_base'], config['fallback_base']):
    if Path(base).is_dir():
        free = shutil.disk_usage(base).free
        candidates.append({'base':base,'free_bytes':free, 'sufficient':free>=config['minimum_free_gib']*1024**3})
selected = next((v['base'] for v in candidates if v['sufficient']), None)
if selected is None:
    raise RuntimeError('No project storage location has sufficient space')
data['storage_candidates'] = candidates
data['remote_work_dir'] = str(Path(selected)/config['project_directory_name'])
data['python_candidates'] = {}
for pat in ('/venv/*/bin/python', '/opt/*/bin/python', '/opt/*/envs/*/bin/python', '/root/*/envs/*/bin/python', '/root/.local/share/mamba/envs/*/bin/python'):
    for py in glob.glob(pat):
        data['python_candidates'][py] = command([py,'-B','-c',"import sys,importlib.util; print(sys.version); print({m: bool(importlib.util.find_spec(m)) for m in ['dolfinx','petsc4py','mpi4py','numpy','scipy','pyvista']})"])
print(json.dumps(data))
"""
            data = json.loads(ssh_python(code))
            write_json(ROOT/"reports"/stage/"remote_environment.json", data)
            (ROOT/"reports"/stage/"remote_environment.txt").write_text(json.dumps(data, indent=2)+"\n")
            evidence["environment"] = data
            evidence["remote_work_dir"] = data["remote_work_dir"]
            print("Remote working directory:", data["remote_work_dir"])
            if a.stage != "0":
                state = json.loads(state_path.read_text())
                if state['remote_work_dir'] != data['remote_work_dir']:
                    raise RuntimeError('Established remote work directory changed; no automatic relocation')
                write_json(ROOT/"reports"/stage/"connection_probe.json", data)
            elif not a.dry_run:
                write_json(state_path, {"host":cfg["host"], "port":cfg["port"],
                           "remote_work_dir": data["remote_work_dir"], "probe_utc":data["timestamp"],
                           "hostname":data["hostname"], "storage_candidates":data["storage_candidates"]})
        else:
            state = json.loads(state_path.read_text())
            if (state["host"], state["port"]) != (cfg["host"], cfg["port"]):
                raise ValueError("Connection changed; run remote_probe.sh again")
            rd = state["remote_work_dir"]
            if rd not in [str(Path(base)/cfg["project_directory_name"]) for base in (cfg["preferred_base"],cfg["fallback_base"])]:
                raise ValueError("Remote directory is outside dedicated project locations")
            evidence["remote_work_dir"] = rd
            print("Remote working directory:", rd, flush=True)
            guard = f"""
from pathlib import Path
p=Path({rd!r})
marker=p/'.fem3d-owner'
if p.is_symlink(): raise RuntimeError('Remote project directory must not be a symlink')
if p.exists() and (not marker.is_file() or marker.read_text().strip()!={OWNER!r}):
    raise RuntimeError('Existing directory is not owned by this project')
"""
            if a.action == "sync":
                code = guard
                if not a.dry_run:
                    code += f"p.mkdir(parents=True,exist_ok=True)\nmarker.write_text({OWNER!r}+'\\n')\n"
                code += "print('Directory ownership check passed')\n"
                print(ssh_python(code))
                manifest = {}
                selected = ["pyproject.toml", "README.md", "src/", "scripts/", "tests/", "configs/", "remote/", "inputs/"]
                for item in selected:
                    path = ROOT/item
                    files = path.rglob("*") if path.is_dir() else [path]
                    for f in files:
                        rel = f.relative_to(ROOT)
                        if f.is_file() and not any(v.startswith(".") or v == "__pycache__" for v in rel.parts) and f.name not in ("connection.local.json", "source_manifest.json"):
                            manifest[str(rel)] = sha256(f)
                write_json(ROOT/"remote/source_manifest.json", {"source_of_truth":str(ROOT), "files":manifest, "timestamp":timestamp(),
                           "wsl_git":command(["git","-C",str(ROOT),"rev-parse","HEAD"])})
                sources = [str(ROOT/v) for v in selected if (ROOT/v).exists()]
                # --relative with /./ keeps the project layout and prevents whole-repo copying.
                sources = [str(ROOT)+"/./"+str(Path(x).relative_to(ROOT)) for x in sources]
                cmd = ["rsync", "-az", "--relative", "--itemize-changes", "--exclude=__pycache__/",
                       "--exclude=.*", "--exclude=connection.local.json", "-e", shlex.join(ssh)]
                if a.dry_run:
                    cmd.append("--dry-run")
                cmd += sources+[cfg["host"]+":"+rd+"/"]
                print(execute(cmd))
            elif a.action == "run":
                if not user_command:
                    p.error("run requires -- command arguments")
                # Validate every synced code/input byte before executing a compute command.
                check = guard + """
import hashlib,json
m=json.loads((p/'remote/source_manifest.json').read_text())
for name,expected in m['files'].items():
    actual=hashlib.sha256((p/name).read_bytes()).hexdigest()
    if actual!=expected: raise RuntimeError('WSL source mismatch: '+name)
print('WSL source hashes verified')
"""
                if a.dry_run:
                    print("Would verify ownership and source hashes before execution")
                else:
                    print(ssh_python(check))
                cmd = [a.python, "-B", "scripts/record_run.py", "--label", a.label,
                       "--stage", str(a.stage), "--ranks", str(a.ranks), "--remote-dir", rd, "--"] + user_command
                remote_command = "cd " + shlex.quote(rd) + " && " + shlex.join(cmd)
                if a.dry_run:
                    print("Would run:", remote_command)
                else:
                    print(execute(ssh+[cfg["host"], remote_command]))
            else:
                if not a.dry_run:
                    print(ssh_python(guard+"print('Directory ownership check passed')\n"))
                dest = ROOT/"outputs"/stage/"remote_return"
                dest.mkdir(parents=True,exist_ok=True)
                cmd = ["rsync", "-az", "--itemize-changes", "--include=/logs/", f"--include=/logs/{stage}/***",
                       "--include=/outputs/", f"--include=/outputs/{stage}/***", "--include=/reports/", f"--include=/reports/{stage}/***",
                       "--exclude=*", "-e", shlex.join(ssh)]
                if a.dry_run:
                    cmd.append("--dry-run")
                cmd += [cfg["host"]+":"+rd+"/", str(dest)+"/"]
                print(execute(cmd))
                evidence["fetched_sha256"] = {str(x.relative_to(dest)):sha256(x) for x in dest.rglob("*") if x.is_file()}
        evidence["returncode"] = 0
    except Exception as exc:
        evidence["error"] = repr(exc)
        evidence["returncode"] = getattr(exc, "returncode", 1)
        raise
    finally:
        evidence["wall_time_s"] = time.perf_counter()-start
        evidence["peak_child_rss_kib"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        evidence["memory_scope"] = "Local SSH/rsync children only; remote compute metadata is fetched separately"
        write_json(log/"metadata.json", evidence)


if __name__ == "__main__":
    main()
