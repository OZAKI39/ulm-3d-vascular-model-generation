from common import *
import subprocess,sys,time,resource
stage=sys.argv[1]
if stage=='prepare':
 while not all((O/'cases'/m/'SUMMARY.json').exists() for m in CASES):time.sleep(15)
 commands=[['analyze.py'],['figures.py'],['render_videos.py','--smoke']]
else:
 assert (O/'SMOKE_VISUAL_REVIEW.json').exists()
 commands=[['render_videos.py'],['report.py'],['final_delivery.py']]
start=time.perf_counter();before=resource.getrusage(resource.RUSAGE_CHILDREN);records=[]
for command in commands:
 t=time.perf_counter();subprocess.run([sys.executable,*command],cwd=A,check=True);records.append(dict(command=command,wall_s=time.perf_counter()-t))
after=resource.getrusage(resource.RUSAGE_CHILDREN)
dump(O/('POSTPROCESS_'+stage.upper()+'_RESOURCES.json'),dict(stage=stage,wall_s=time.perf_counter()-start,child_cpu_s=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,peak_child_RSS_bytes=after.ru_maxrss*1024,commands=records))
print('POSTPROCESS_STAGE_DONE',stage,flush=True)
