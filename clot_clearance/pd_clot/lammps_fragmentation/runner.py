"""Run one real case with MPI and measured process-tree resources."""
from pathlib import Path
import os,subprocess,time,re
import psutil
from .config import LMP,dump

def run_lammps(case):
    p=Path(case['path']);r=case['request']
    if (p/'result.json').exists():raise FileExistsError('Successful case is immutable; choose a new case directory')
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    cmd=['mpirun','--bind-to','core','-np',str(r['mpi_ranks']),str(LMP),'-in','in.lammps','-log','log.lammps','-screen','screen.lammps']
    command=['/usr/bin/time','-f','%U %S %e %M %x','-o','resource.time',*cmd];dump(p/'RUN_COMMAND.json',dict(argv=cmd,cwd=p,environment={k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']}))
    start=time.perf_counter();peak=0;samples=[]
    with (p/'launcher.log').open('w') as log:
        proc=subprocess.Popen(command,cwd=p,env=env,stdout=log,stderr=subprocess.STDOUT);parent=psutil.Process(proc.pid)
        while proc.poll() is None:
            try:
                children=[parent,*parent.children(recursive=True)];rss=sum(q.memory_info().rss for q in children if q.is_running());peak=max(peak,rss)
                samples.append((time.perf_counter()-start,rss))
            except (psutil.NoSuchProcess,psutil.AccessDenied):pass
            time.sleep(.2)
    wall=time.perf_counter()-start;numbers=(p/'resource.time').read_text().splitlines()[-1].split();ucpu,scpu,elapsed,rss,code=map(float,numbers)
    output=(p/'log.lammps').read_text() if (p/'log.lammps').exists() else '';loop=re.search(r'Loop time of ([\d.e+-]+) on (\d+) procs for (\d+) steps',output)
    timings={}
    for line in output.splitlines():
        fields=[a.strip() for a in line.split('|')]
        if len(fields)>=5 and fields[0] in ['Pair','Bond','Neigh','Comm','Output','Modify','Other']:
            try:timings[fields[0]]=dict(min_s=float(fields[1]),avg_s=float(fields[2]),max_s=float(fields[3]),percent=float(fields[-1]))
            except ValueError:pass
    perf=dict(exit_code=proc.returncode,wall_seconds=wall,CPU_seconds=ucpu+scpu,peak_RSS_bytes=max(peak,int(rss)*1024),
        mpi_ranks=r['mpi_ranks'],node_count=case['geometry']['info']['node_count'],initial_bond_count=case['geometry']['info']['bond_count'],
        average_bonds_per_node=case['geometry']['info']['average_bonds_per_node'],explicit_steps=int(loop[3]) if loop else 0,
        loop_wall_s=float(loop[1]) if loop else None,steps_per_second=int(loop[3])/float(loop[1]) if loop else None,
        LAMMPS_sections=timings,memory_accounting='Maximum sampled process-tree RSS and GNU time individual high-water mark',scope='Actual hemisphere production case')
    dump(p/'PERFORMANCE.json',perf)
    with (p/'resource_samples.csv').open('w') as f:
        f.write('elapsed_s,process_tree_RSS_bytes\n');f.writelines(f'{t},{m}\n' for t,m in samples)
    if proc.returncode or 'ERROR:' in output:raise RuntimeError(f'LAMMPS failed; inspect {p}/log.lammps and launcher.log')
    return perf
