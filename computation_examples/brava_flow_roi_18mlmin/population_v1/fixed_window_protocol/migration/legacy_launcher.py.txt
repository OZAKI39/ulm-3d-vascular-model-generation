"""New campaign only; records process identity and GPU/CPU telemetry."""
from pathlib import Path
import argparse,csv,json,os,subprocess,time,sys
import psutil
ap=argparse.ArgumentParser();ap.add_argument('--git-commit',required=True);ap.add_argument('--run-name',default='realization_0000');args=ap.parse_args()
root=Path(__file__).resolve().parents[1];output=root/'runs'/args.run_name
if output.exists():raise SystemExit('Refusing to overwrite existing production output')
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
    PYTHONPATH=str(root/'source/particle_3d/src'))
command=[str(root/'env/bin/python'),'-u','-B','-m','particle_3d.population_v1.run',
    '--environment',str(root/'data/environment.json'),'--output',str(output),
    '--backend','cuda','--concentration','4e9','--seed','2026092909','--minimum-exits','20',
    '--post-infusion','12','--git-commit',args.git_commit]
log=(root/'logs'/f'{args.run_name}.log').open('w')
process=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT)
receipt=dict(wrapper_pid=os.getpid(),solver_pid=process.pid,command=command,started_unix_s=time.time(),
    scope='NEW_BRAVA_GLOBAL_TIME_POPULATION_ONLY_NO_CFD_NO_OLD_STOPPED_BATCH')
(root/'logs'/f'{args.run_name}_PROCESS.json').write_text(json.dumps(receipt,indent=2)+'\n')
proc=psutil.Process(process.pid);proc.cpu_percent()
with (root/'logs'/f'{args.run_name}_telemetry.csv').open('w',newline='') as stream:
    writer=csv.writer(stream);writer.writerow(['unix_s','gpu_utilization_percent','gpu_memory_MiB','gpu_power_W','process_CPU_percent_one_core_100','process_RSS_bytes'])
    while process.poll() is None:
        try:
            values=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,power.draw','--format=csv,noheader,nounits'],text=True,timeout=5).strip().split(',')
            writer.writerow([time.time(),*[v.strip() for v in values],proc.cpu_percent(),proc.memory_info().rss]);stream.flush()
        except (subprocess.SubprocessError,psutil.Error):pass
        time.sleep(.5)
receipt.update(exit_code=process.wait(),ended_unix_s=time.time());log.close()
(root/'logs'/f'{args.run_name}_PROCESS.json').write_text(json.dumps(receipt,indent=2)+'\n')
raise SystemExit(receipt['exit_code'])
