#!/usr/bin/env python3
"""Userspace process/cgroup telemetry; never captures credentials or env vars."""
from pathlib import Path
import argparse,json,os,socket,time
import psutil


def sample():
    rows=[]
    for p in psutil.process_iter(['pid','name','cmdline','memory_info','cpu_times']):
        try:
            command=' '.join(p.info['cmdline'] or [])
            if 'particle82a_' not in command or p.pid==os.getpid() or p.info['name'] not in ['python','python3']:
                continue
            cpu=p.info['cpu_times'];m=p.memory_full_info()
            rows.append(dict(pid=p.pid,cpu_seconds=cpu.user+cpu.system,rss_bytes=m.rss,pss_bytes=m.pss))
        except (psutil.Error,TypeError):pass
    cgroup={}
    for name in ['memory.current','memory.max','cpu.stat','cpu.max']:
        p=Path('/sys/fs/cgroup')/name
        if p.exists():cgroup[name]=p.read_text().strip()
    return dict(time=time.time(),hostname=socket.gethostname(),processes=rows,
        total_process_pss_bytes=sum(p['pss_bytes'] for p in rows),cgroup=cgroup,
        system_available_memory_bytes=psutil.virtual_memory().available)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--stop-file',required=True)
    a=p.parse_args();out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('a') as f:
        while not Path(a.stop_file).exists():
            f.write(json.dumps(sample())+'\n');f.flush();time.sleep(10)


if __name__=='__main__':main()
