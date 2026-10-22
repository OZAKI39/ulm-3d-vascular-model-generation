from pathlib import Path
import sys,json,subprocess,hashlib,importlib.metadata as m,time
p=Path('/workspace/mevo_acoustic_microbubbles/f10z_simulation_20261018/results/straight_tube_exp1/f10z/evidence');p.mkdir(exist_ok=True)
versions={}
for name in ['numpy','scipy','h5py','osqp','warp-lang','PyYAML','psutil','torch']:
 try:versions[name]=m.version(name)
 except m.PackageNotFoundError:versions[name]='NOT_INSTALLED'
r={'time_unix':time.time(),'python':sys.version,'prefix':sys.prefix,'dependencies':versions}
for name,args in [('gpu',['nvidia-smi','--query-gpu=name,memory.total,memory.used,utilization.gpu,driver_version','--format=csv']),('cpu',['lscpu']),('memory',['free','-h']),('disk',['df','-h','/workspace']),('active_processes',['ps','-eo','pid,ppid,pcpu,rss,args'])]:
 out=subprocess.check_output(args,text=True)
 if name=='active_processes':out='\n'.join(line for line in out.splitlines() if 'f10z.run' in line or 'clot_biphasic' in line)
 r[name]=out
(p/'environment_remote.json').write_text(json.dumps(r,indent=2));(p/'environment_remote.lock.txt').write_text('\n'.join(sorted(d.metadata['Name']+'=='+d.version for d in m.distributions()))+'\n');print(json.dumps({'dependencies':versions,'gpu':r['gpu'],'f10z_processes':r['active_processes']}))
