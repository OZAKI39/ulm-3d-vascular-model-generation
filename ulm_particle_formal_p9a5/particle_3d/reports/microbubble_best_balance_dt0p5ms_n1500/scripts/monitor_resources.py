"""One-second observations; stop with the campaign completion file."""
from pathlib import Path
import csv,json,subprocess,time
import psutil
HERE=Path(__file__).resolve().parents[1]

with (HERE/'logs/resource_usage.csv').open('w',newline='') as stream:
    fields=['unix_s','cpu_percent_host_observation','process_cpu_s','process_rss_bytes',
            'gpu_utilization_percent','gpu_memory_used_MiB','gpu_power_W']
    w=csv.DictWriter(stream,fieldnames=fields);w.writeheader()
    while not (HERE/'data/final_summary.json').exists():
        rss=cpu=0
        for p in psutil.process_iter(['cmdline','memory_info','cpu_times']):
            try:
                if str(HERE) not in ' '.join(p.info['cmdline'] or []):continue
                rss+=p.info['memory_info'].rss;cpu+=sum(p.info['cpu_times'][:2])
            except (psutil.NoSuchProcess,psutil.AccessDenied):pass
        gpu=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,power.draw',
                                     '--format=csv,noheader,nounits'],text=True).strip().split(',')
        w.writerow(dict(zip(fields,[time.time(),psutil.cpu_percent(),cpu,rss,*map(float,gpu)])))
        stream.flush();time.sleep(1)
