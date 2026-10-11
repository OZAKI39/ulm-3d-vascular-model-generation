import os,subprocess,time
from pathlib import Path

def sample(root_pid):
 info={}
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   raw=(p/'stat').read_text();v=raw[raw.rfind(')')+2:].split();info[int(p.name)]=(int(v[1]),int(v[21])*os.sysconf('SC_PAGE_SIZE'),(int(v[11])+int(v[12]))/os.sysconf('SC_CLK_TCK'))
  except (OSError,ValueError,IndexError):pass
 ids={root_pid}
 for _ in range(8):
  old=len(ids);ids.update(k for k,v in info.items() if v[0] in ids)
  if len(ids)==old:break
 gpu=None
 try:gpu=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'],text=True,timeout=3).strip()
 except Exception:pass
 return dict(unix_time=time.time(),processes=len(ids),RSS_bytes=sum(info[p][1] for p in ids if p in info),CPU_seconds_live=sum(info[p][2] for p in ids if p in info),gpu_util_percent_memory_MiB=gpu)
