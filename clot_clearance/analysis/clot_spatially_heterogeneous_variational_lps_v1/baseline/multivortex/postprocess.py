from common import *
import time,subprocess,sys,resource
started=time.perf_counter()
while not all((O/'cases'/m/'SUMMARY.json').exists() for m in CASES):time.sleep(10)
for script in ['analyze.py','figures.py','report.py']:
    print('POSTPROCESS_START',script,flush=True)
    subprocess.run([sys.executable,str(A/script)],check=True,cwd=R)
usage=resource.getrusage(resource.RUSAGE_CHILDREN)
dump(O/'POSTPROCESS_RESOURCES.json',dict(wall_including_wait_s=time.perf_counter()-started,child_cpu_s=usage.ru_utime+usage.ru_stime,child_peak_RSS_bytes=usage.ru_maxrss*1024,waiting_is_not_computing=True))
print('POSTPROCESS_DONE',flush=True)
