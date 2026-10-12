"""Read-only final environment and supervisor process timing capture."""

from datetime import datetime
import json
import re
import shlex
import subprocess
from .common import BASE, write


def main():
    code = """import json,sys,platform,importlib.metadata,subprocess
from pathlib import Path
print(json.dumps(dict(python=sys.version,cpu=next(x.split(':',1)[1].strip() for x in Path('/proc/cpuinfo').read_text().splitlines() if x.startswith('model name')),platform=platform.platform(),packages={k:importlib.metadata.version(k) for k in ['warp-lang','numpy','scipy','h5py']},gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total,utilization.gpu','--format=csv,noheader'],text=True).strip(),supervisor=subprocess.run(['supervisorctl','status'],text=True,capture_output=True).stdout,log=[line for line in Path('/var/log/supervisor/supervisord.log').read_text().splitlines() if 'mevo_fast_pzt_f9_' in line])))"""
    data = json.loads(
        subprocess.check_output(
            [
                "ssh",
                "mevo-new",
                "/workspace/mevo_acoustic_microbubbles/.venv/bin/python -c " + shlex.quote(code),
            ],
            text=True,
        )
    )
    log = data.pop("log")
    (BASE / "SUPERVISOR_RECORD.log").write_text("\n".join(log) + "\n")
    jobs = {}
    for line in log:
        label = re.search(r"mevo_fast_pzt_f9_(\w+)", line)[1]
        time = datetime.strptime(line[:23], "%Y-%m-%d %H:%M:%S,%f")
        if "spawned:" in line:
            jobs.setdefault(label, {})["start"] = time.isoformat()
        if "exited:" in line:
            row = jobs[label]
            assert "exit status 0; expected" in line
            row.update(
                stop=time.isoformat(),
                exit_record=line,
                process_elapsed_s=(time - datetime.fromisoformat(row["start"])).total_seconds(),
            )
    assert len(jobs) == 6 and all("stop" in v for v in jobs.values())
    assert not any(
        "mevo_fast_pzt_f9_" in line and "RUNNING" in line for line in data["supervisor"].splitlines()
    )
    write(
        BASE / "JOB_TIMING.json",
        dict(
            jobs=jobs,
            total_process_elapsed_s=sum(v["process_elapsed_s"] for v in jobs.values()),
            scope="Supervisor process start-to-exit sums for six serial F9 jobs. Includes setup, warmup, CPU/GPU and output; excludes development, WSL analysis/plots and transfers. Not exclusive GPU time.",
        ),
    )
    write(
        BASE / "ENVIRONMENT.json",
        dict(
            **data,
            no_dependency_upgrades=True,
            timestamp=datetime.now().astimezone().isoformat(),
            server_root="/workspace/mevo_acoustic_microbubbles/f9",
        ),
    )


if __name__ == "__main__":
    main()
