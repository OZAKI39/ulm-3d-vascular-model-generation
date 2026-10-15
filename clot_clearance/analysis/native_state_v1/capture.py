from pathlib import Path
import json,hashlib,subprocess,sys,os,shutil
C=json.loads(Path(__file__).with_name('config.json').read_text());ROOT=Path.cwd();OUT=ROOT/'verification'/C['task']/C['run'];P=OUT/'provenance';P.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cmd(s):return subprocess.run(s,shell=True,text=True,capture_output=True).stdout
src=ROOT/'tools/lammps-contact-onset-v1';build=ROOT/'tools/lammps-native-y-wall-shared-v1-build'
files=['src/BPM/bond_bpm_peri.cpp','src/BPM/bond_bpm_peri.h','src/BPM/bond_bpm.cpp','src/BPM/bond_bpm.h','src/BPM/pair_bpm_peri.cpp','src/atom.cpp','src/atom_vec.cpp','src/comm_brick.cpp','src/fix_property_atom.cpp','src/fix_bond_history.cpp','src/verlet.cpp']
ids={}
for f in files:
 p=src/f;dest=P/'native_source'/f;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest);ids[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
for f in ['liblammps.so','lmp','CMakeCache.txt','CMakeFiles/lammps.dir/flags.make']:
 p=build/f
 if p.exists():ids[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
 if p.suffix in ['.txt','.make']:shutil.copyfile(p,P/p.name)
for f in ['AGENTS.md','PROJECT_ACTIVE_CONTEXT.md']:
 if (ROOT/f).exists():shutil.copyfile(ROOT/f,P/f)
identity={'files':ids,'pwd':str(ROOT),'git_status':cmd('git status --short'),'branch':cmd('git branch --show-current'),'head':cmd('git rev-parse HEAD'),'native_git_head':cmd(f'git -C {src} rev-parse HEAD'),'native_git_status':cmd(f'git -C {src} status --short'),'python':sys.version,'python_path':sys.executable,'mpi':cmd('which mpiexec; mpiexec --version'),'compiler':cmd('g++ --version'),'cpu':cmd('lscpu'),'cgroup':cmd('cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us /sys/fs/cgroup/cpu/cpu.cfs_period_us /sys/fs/cgroup/memory/memory.limit_in_bytes'),'memory':cmd('free -b'),'jobs':cmd('ps -eo pid,comm,pcpu,pmem --sort=-pcpu | head -25'),'disk':cmd('df -B1 /workspace')}
import lammps
l=lammps.lammps(cmdargs=['-log','none','-screen','none']);identity.update(lammps_version=l.version(),lammps_python=lammps.__file__,loaded_libraries=[x.strip() for x in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in x or 'libmpi.so' in x]);l.close()
(OUT/'SOURCE_AND_BINARY_IDENTITY.json').write_text(json.dumps(identity,indent=2))
# All existing verification/results scientific files, excluding this task. Hash also exact source and production binary.
protected={}
for folder in ['verification','results']:
 for p in (ROOT/folder).rglob('*'):
  if p.is_file() and C['task'] not in p.parts:protected[str(p.relative_to(ROOT))]=sha(p)
protected.update({str(Path(p).relative_to(ROOT)):v['sha256'] for p,v in ids.items()})
(P/'PROTECTED_START.json').write_text(json.dumps(protected,indent=2))
print(json.dumps({'version':identity['lammps_version'],'protected_files':len(protected),'cgroup':identity['cgroup'],'jobs':identity['jobs'],'memory':identity['memory']}))
