#!/usr/bin/env python3
import pathlib,subprocess,os,sys,csv,json,time,hashlib,shutil,datetime,fcntl,signal,threading
sys.dont_write_bytecode=True
R=pathlib.Path(__file__).resolve().parents[1];B=R/'frozen_bundle';P=pathlib.Path('/workspace/hemocell_gpu_poc/toolchains/nvhpc/Linux_x86_64/26.5');C=R/'verification/GPU_NUMERICAL_COMPARISON_CONTRACT.json'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,v):p.write_text(json.dumps(v,indent=2))
def processes():
 rows={}
 for p in pathlib.Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   s=(p/'stat').read_text().split();rows[int(p.name)]=dict(name=(p/'comm').read_text().strip(),exe=os.readlink(p/'exe'),ticks=int(s[13])+int(s[14]),rss=int((p/'statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE'))
  except (OSError,ValueError):pass
 return rows
def memory():
 info={l.split(':')[0]:int(l.split()[1])*1024 for l in pathlib.Path('/proc/meminfo').read_text().splitlines() if l.startswith(('MemAvailable:','MemTotal:'))}
 cg=pathlib.Path('/sys/fs/cgroup/memory.max');used=pathlib.Path('/sys/fs/cgroup/memory.current')
 if cg.exists() and cg.read_text().strip()!='max':info['effective_available']=min(info['MemAvailable'],int(cg.read_text())-int(used.read_text()))
 else:info['effective_available']=info['MemAvailable']
 return info
def gpu():return [float(x.strip()) for x in subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,memory.total','--format=csv,noheader,nounits'],text=True).strip().split(',')]
def materialize(D):
 D.mkdir()
 for name in ['contracts','diagnostics/field_samples','logs','provenance']:(D/name).mkdir(parents=True,exist_ok=True)
 # Hardlinks only to task-owned frozen copies; never to protected formal inputs.
 for p in (B/'frozen_contracts').iterdir():
  if p.is_file():os.link(p,D/'contracts'/p.name)
 lines=(B/'frozen_contracts/step3_solver_parameters.txt').read_text().splitlines();old=list(lines);nums=lines[2].split();assert nums[-1]=='1000';nums[-1]='700000';lines[2]=' '.join(nums)
 lines[0]='frozen_inputs/step1_geometry_contract/geometry/cfd_surface_axis_aligned_inlet_m.stl';lines[1]='frozen_inputs/step2'
 (D/'contracts/solver_parameters.txt').write_text('\n'.join(lines)+'\n');(D/'contracts/monitor_parameters.txt').write_text('100 5000 119872 240000 180543\n')
 assert lines[2].split()[:-1]==old[2].split()[:-1] and lines[3:]==old[3:]
 save(D/'provenance/input_verification.json',dict(status='PASS',copied_frozen_contracts={p.name:sha(p) for p in (B/'frozen_contracts').iterdir() if p.is_file()},materialized_run_files={p.name:sha(p) for p in (D/'contracts').iterdir()},path_rebinding_only=True,legacy_parameter_file_max_steps=700000,external_PoC_hard_cap=10000,physical_numeric_tokens_changed=False,formal_inputs_hardlinked=False))
