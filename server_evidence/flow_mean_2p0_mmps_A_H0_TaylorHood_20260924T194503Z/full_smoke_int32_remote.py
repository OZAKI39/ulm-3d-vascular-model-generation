"""One-step GPU smoke; edits only fresh per-attempt configuration files."""
from pathlib import Path
import argparse,json,subprocess,xml.etree.ElementTree as ET,os,time,re,signal
p=argparse.ArgumentParser();p.add_argument('root');p.add_argument('--case',default='case_int32');p.add_argument('--shift',action='store_true');a=p.parse_args()
root=Path(a.root);b=json.loads((root/'audit/remote_protection_manifest.json').read_text())['build']
name='smoke_shift' if a.shift else 'smoke_baseline';out=root/a.case/name;out.mkdir()
tree=ET.parse(root/a.case/'run/solver.xml')
changes={'.//Number_of_time_steps':'1','.//Increment_in_saving_VTK_files':'1',
         './/Increment_in_saving_restart_files':'1','.//Add_equation/Min_iterations':'1',
         './/Add_equation/Max_iterations':'1'}
for key,value in changes.items():tree.find(key).text=value
tree.write(out/'solver.xml',encoding='utf-8',xml_declaration=True)
opts=(root/'baseline_PETSC_OPTIONS.txt').read_text().strip()
extra=' -sub_pc_factor_shift_type nonzero -sub_pc_factor_shift_amount 1e-6' if a.shift else ''
opts+=extra;(out/'PETSC_OPTIONS.txt').write_text(opts+'\n')
(out/'config_diff.json').write_text(json.dumps(dict(smoke_duration_only=changes,solver_only_extra_options=extra),indent=2)+'\n')
# ABI offsets were independently compiled against the unchanged actual build headers.
gdb='''set pagination off
set confirm off
set breakpoint pending on
break *'fs::get_thood_fs(ComMod&, std::array<fsType, 2ul>&, mshType const&, bool, int)'
commands
silent
set $mesh = (char*)$rdx
set $fsp = *(char**)($mesh+896)
printf "TH_RUNTIME mesh_eNoN=%d nFs=%d mesh_eType=%d velocity_eNoN=%d pressure_eNoN=%d velocity_eType=%d pressure_eType=%d vmsStab=%d option=%d\\n", *(int*)($mesh+40), *(int*)($mesh+64), *(int*)($mesh+36), *(int*)($fsp+8), *(int*)($fsp+208), *(int*)($fsp+4), *(int*)($fsp+204), (int)$rcx, (int)$r8
disable 1
continue
end
run
'''
(out/'runtime_audit.gdb').write_text(gdb)
env={k:os.environ[k] for k in ('HOME','USER','LOGNAME','LANG') if k in os.environ}
env.update(PATH='/usr/bin:/bin:/usr/local/cuda/bin',LC_ALL='C',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',LD_LIBRARY_PATH=b['runtime_library_path'],PETSC_OPTIONS=opts)
w=b['PETSc_build'];cmd=[w['candidate_wrapper'],w['launcher'],'-n','1',w['candidate_wrapper'],'/usr/bin/gdb','--batch','-x','runtime_audit.gdb','--args',b['executable'],'solver.xml']
start=time.time();samples=[];stop_reason=None
with (out/'solver.log').open('w') as f:
    proc=subprocess.Popen(cmd,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
    while proc.poll() is None:
        sample={'elapsed_s':time.time()-start}
        smi=subprocess.run(['nvidia-smi','--query-gpu=memory.used,utilization.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True)
        sample['nvidia_smi']=smi.stdout.strip()
        host=[]
        for d in Path('/proc').iterdir():
            if not d.name.isdigit():continue
            try:
                if str(out) not in os.readlink(d/'cwd'):continue
                status=(d/'status').read_text();name=re.search(r'^Name:\s*(.*)$',status,re.M).group(1)
                rss=re.search(r'^VmRSS:\s*(\d+)',status,re.M);hwm=re.search(r'^VmHWM:\s*(\d+)',status,re.M)
                if rss:host.append(dict(pid=int(d.name),name=name,rss_kib=int(rss.group(1)),hwm_kib=int(hwm.group(1))))
            except (OSError,PermissionError,AttributeError):pass
        sample['processes']=host;samples.append(sample)
        if time.time()-start>1800:
            stop_reason='smoke_wall_time_limit';os.killpg(proc.pid,signal.SIGTERM);break
        (out/'resource_samples.json').write_text(json.dumps(samples)+'\n');time.sleep(3)
    try:proc.wait(timeout=30)
    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
log=(out/'solver.log').read_text(errors='replace')
record=dict(command=cmd,exit_code=proc.returncode,elapsed_seconds=time.time()-start,stop_reason=stop_reason,
  solver_only_extra_options=extra,runtime_lines=[l for l in log.splitlines() if 'TH_RUNTIME mesh_' in l],
  linear_solves=len(re.findall(r'SV13Q_END',log)),linear_failures=len(re.findall(r'SV13Q_END[^\n]*reason=-',log)),
  gpu_matrix='seqaijcusparse' in log,gpu_vector='seqcuda' in log,
  results=[str(x) for x in out.rglob('result*.vtu')],samples=samples)
(out/'execution.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k!='samples'},indent=2),flush=True)
print(log[-3500:],flush=True)
