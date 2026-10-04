import psutil,time,json
from pathlib import Path
out=Path(__file__).resolve().parent.parent;rows=[];peak=0;least=None;swapped=0
while not (out/'provenance/STOP_RESOURCE_MONITOR').exists():
 roots=[]
 for p in psutil.process_iter(['pid','cmdline']):
  try:
   args=p.info['cmdline'] or []
   if '-m' in args and any(v in args for v in ['scripts.c7_early_branch','scripts.c7_spatial_evidence']):roots.append(p)
  except psutil.Error:pass
 processes={p.pid:p for root in roots for p in [root]+root.children(recursive=True)};rss=0;swap=0
 for pid,p in processes.items():
  try:
   rss+=p.memory_info().rss
   for s in Path(f'/proc/{pid}/status').read_text().splitlines():
    if s.startswith('VmSwap:'):swap+=int(s.split()[1])*1024
  except (psutil.Error,FileNotFoundError,ProcessLookupError):pass
 mem=psutil.virtual_memory()
 if roots:
  peak=max(peak,rss);least=mem.available if least is None else min(least,mem.available);swapped=max(swapped,swap)
  rows.append(dict(time=time.time(),processes=len(processes),tree_RSS_bytes=rss,tree_swap_bytes=swap,RAM_available=mem.available))
  result=dict(peak_observed_tree_RSS_bytes=peak,min_observed_RAM_available=least,total_RAM=mem.total,max_observed_task_swap_bytes=swapped,samples=rows,scope='External synchronized samples, started during initial C7 broad scan; prior C5 stages report separate per-process peaks')
  temp=out/'provenance/RESOURCE_MONITOR.tmp';temp.write_text(json.dumps(result,indent=2));temp.replace(out/'provenance/RESOURCE_MONITOR.json')
 time.sleep(1)
