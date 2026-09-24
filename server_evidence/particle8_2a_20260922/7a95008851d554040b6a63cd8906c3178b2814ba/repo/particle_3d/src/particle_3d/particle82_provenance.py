"""Explicit host receipts for P8.2. No credentials or full environment capture."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import platform
import shutil
import socket
import subprocess
import sys


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    os.replace(temp, path)


def command(args):
    try:
        result = subprocess.run(args, text=True, capture_output=True, timeout=45)
        return dict(command=args, exit_code=result.returncode, stdout=result.stdout, stderr=result.stderr)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return dict(command=args, exit_code=None, error=str(exc))


def collect_host(role):
    commands = [['hostname'], ['uname', '-a'], ['nproc'], ['lscpu'], ['free', '-h'], ['df', '-h'],
                ['python', '--version'], ['which', 'python'], ['python3', '--version'], ['which', 'python3'],
                ['git', 'rev-parse', 'HEAD'], ['bash', '-c', 'ulimit -a']]
    if shutil.which('nvidia-smi'):
        commands += [['nvidia-smi'], ['nvidia-smi', '--query-gpu=name,memory.total,driver_version', '--format=csv,noheader']]
    mem = {line.split(':')[0]: line.split(':')[1].strip() for line in Path('/proc/meminfo').read_text().splitlines()}
    limits = {}
    for name in ['memory.max', 'memory.current', 'cpu.max', 'cpuset.cpus.effective']:
        p = Path('/sys/fs/cgroup') / name
        if p.exists(): limits[name] = p.read_text().strip()
    return dict(role=role, captured_utc=datetime.now(timezone.utc).isoformat(), hostname=socket.gethostname(),
                python=sys.version, executable=sys.executable, platform=platform.platform(),
                nproc=len(os.sched_getaffinity(0)), os_cpu_count=os.cpu_count(), memory=mem, cgroup=limits,
                gpu_solver_status='GPU_AVAILABLE_BUT_NOT_USED_BY_CURRENT_SOLVER' if shutil.which('nvidia-smi') else 'GPU_NOT_AVAILABLE',
                commands=[command(args) for args in commands])


def require_remote(receipt, expected_hostname):
    if receipt.get('role') != 'HEAVY_COMPUTE_PRIMARY':
        raise ValueError('Formal trajectory computation requires HEAVY_COMPUTE_PRIMARY')
    if receipt.get('hostname') != expected_hostname or socket.gethostname() != expected_hostname:
        raise ValueError('Actual compute hostname does not match SSH-attested remote hostname')
    if 'microsoft' in platform.release().lower() or 'wsl' in platform.release().lower():
        raise ValueError('Formal heavy trajectory computation is forbidden on WSL')
    if receipt.get('ssh_endpoint') != 'root@50.115.148.16:4159':
        raise ValueError('Missing authorized remote endpoint attestation')
    return True


if __name__ == '__main__':
    print(json.dumps(collect_host(sys.argv[1]), ensure_ascii=False))
