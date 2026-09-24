#!/usr/bin/env python3
"""Capture local/SSH host provenance before any new heavy computation."""
from pathlib import Path
import json
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from particle_3d.particle82_provenance import collect_host, atomic_json, sha256

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / 'particle_3d/reports/particle8_2'
SSH = ['ssh', '-p', '4159', '-i', str(Path.home()/'.ssh/vast_step3b_ed25519'),
       '-o', 'IdentitiesOnly=yes', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15',
       '-o', 'ServerAliveInterval=30', 'root@50.115.148.16']


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    local = collect_host('LOCAL_WSL_DEVELOPMENT_HOST')
    local['compute_policy'] = 'DEVELOPMENT_ORCHESTRATION_ONLY; at most 10 trajectory smoke; excluded from formal catalog'
    atomic_json(ROOT/'LOCAL_HOST_PROVENANCE.json', local)
    source = (REPO/'particle_3d/src/particle_3d/particle82_provenance.py').read_text()
    result = subprocess.run(SSH+['python3 - HEAVY_COMPUTE_PRIMARY'], input=source, text=True, capture_output=True, timeout=120)
    (ROOT/'logs').mkdir(exist_ok=True)
    (ROOT/'logs/ssh_provenance_stderr.log').write_text(result.stderr)
    if result.returncode:
        atomic_json(ROOT/'REMOTE_CONNECTION_FAILURE.json', dict(exit_code=result.returncode, stderr=result.stderr))
        raise RuntimeError('Remote authentication/probe failed; no heavy local fallback permitted')
    remote = json.loads(result.stdout)
    remote['ssh_endpoint'] = 'root@50.115.148.16:4159'
    remote['attestation_method'] = 'COMMAND_EXECUTED_THROUGH_AUTHENTICATED_SSH_FROM_LOCAL_ORCHESTRATOR'
    atomic_json(ROOT/'REMOTE_HOST_PROVENANCE.json', remote)
    old = REPO/'particle_3d/outputs/particle8_1'
    inspected = [old/'data/software_environment.json', old/'data/test_run.json', old/'data/run_summary.json',
                 old/'ARTIFACT_SHA256.json', old/'OUTPUT_FILES.json']
    inspected += sorted((old/'logs').glob('*.log')) + sorted((old/'trajectories').glob('*.json'))
    host_keys = {'hostname', 'compute_host', 'server_hostname', 'server_ip', 'worker_hostname'}
    found = []
    def scan(value, filename, context=''):
        if isinstance(value, dict):
            for key, val in value.items():
                if key.lower() in host_keys: found.append(dict(path=filename, key=context+key, value=val))
                scan(val, filename, context+key+'.')
        elif isinstance(value, list):
            for val in value: scan(val, filename, context)
    for p in inspected:
        if p.suffix == '.json': scan(json.loads(p.read_text()), str(p.relative_to(REPO)))
    prior = dict(PARTICLE8_1_HEAVY_COMPUTE_HOST='NOT_PROVEN',
        reason='No authenticated per-trajectory/shard compute-host receipt. The recorded WSL platform may describe export/render/test environment; it does not prove where every heavy integration ran.',
        environment_platform=json.loads((old/'data/software_environment.json').read_text())['platform'],
        workers_in_run_summary=json.loads((old/'data/run_summary.json').read_text())['workers'],
        explicit_host_fields=found, inspected_sha256={str(p.relative_to(REPO)):sha256(p) for p in inspected})
    doc=dict(stage='Particle-8.2', baseline_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=REPO).strip(),
        local=local, remote=remote, particle8_1=prior,
        policy=dict(local='DEVELOPMENT_ORCHESTRATION_ONLY',remote='HEAVY_COMPUTE_PRIMARY',local_smoke_limit=10))
    atomic_json(ROOT/'COMPUTE_PROVENANCE.json',doc)
    (ROOT/'COMPUTE_PROVENANCE.md').write_text('\n'.join([
        '# Particle-8.2 计算来源审计', '',
        f'本地：{local["hostname"]}，可用逻辑 CPU={local["nproc"]}，RAM={local["memory"]["MemTotal"]}。职责仅开发、轻量测试、至多 10 条 smoke、调度和审核。',
        f'远程：{remote["hostname"]}，SSH root@50.115.148.16:4159，可用逻辑 CPU={remote["nproc"]}，RAM={remote["memory"]["MemTotal"]}。正式大计算必须在该主机执行；容器内存/CPU 上限另见 JSON。',
        '', '**PARTICLE8_1_HEAVY_COMPUTE_HOST = NOT_PROVEN**。P8.1 环境文件记录了 WSL 平台，但没有逐轨迹/分片的主机身份与 SSH 计算回执。这不足以证明正式重积分在远程执行，也不能仅据此断言全部重积分发生在 WSL。',
        f'已检查 {len(inspected)} 个历史元数据、日志和清单，保存其 SHA256；历史文件不改写。',
        '', '完整 hostname、uname、nproc、lscpu、free、df、Python、git、ulimit 和 GPU 查询原始输出见 COMPUTE_PROVENANCE.json。',
        '后续正式分片须携带主机、PID、worker、源码提交、配置与 Frozen 输入散列；任何本地正式分片都会使验证失败。', '']))
    print(json.dumps(dict(local_hostname=local['hostname'],remote_hostname=remote['hostname'],local_nproc=local['nproc'],remote_nproc=remote['nproc'],remote_ram=remote['memory']['MemTotal'],prior='NOT_PROVEN')))


if __name__=='__main__':main()
