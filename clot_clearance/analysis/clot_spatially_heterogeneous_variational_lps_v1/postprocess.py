from common import *
import subprocess,sys,time,resource,importlib.metadata,argparse
stages=['literature.py','analyze.py','initial_history_audit.py','figures.py','render_videos.py','report.py','final_delivery.py']
parser=argparse.ArgumentParser();parser.add_argument('--resume-from',choices=stages,default=stages[0]);args=parser.parse_args()
start=time.perf_counter();cpu=time.process_time()
while not all((O/'cases'/m/'SUMMARY.json').exists() for m in CASES):
 print('WAITING_FOR_FORMAL_NATIVE_COMPLETION',flush=True);time.sleep(20)
assert all(json.loads((O/'cases'/m/'SUMMARY.json').read_text())['status']=='COMPLETED' for m in CASES)
start=time.perf_counter();cpu=time.process_time()
packages={}
for k in ['numpy','scipy','gstools','gstools-cython','pandas','matplotlib','pyvista','vtk','mpi4py','pillow']:
 try:packages[k]=importlib.metadata.version(k)
 except importlib.metadata.PackageNotFoundError:packages[k]=None
dump(O/'ENVIRONMENT.json',dict(python=sys.version,packages=packages,threads=dict(OMP_NUM_THREADS=os.getenv('OMP_NUM_THREADS'),OPENBLAS_NUM_THREADS=os.getenv('OPENBLAS_NUM_THREADS')),ffmpeg=subprocess.check_output(['ffmpeg','-version'],text=True).splitlines()[0]))
import shutil
shutil.copy2(A/'REPRODUCE.md',O/'REPRODUCE.md')
for script in stages[stages.index(args.resume_from):]:
 print('POSTPROCESS',script,flush=True);subprocess.run([sys.executable,'-u',str(A/script)],check=True)
dump(O/'POSTPROCESS_RESOURCES.json',dict(wall_s=time.perf_counter()-start,driver_cpu_s=time.process_time()-cpu,children_cpu_s=resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime,peak_child_RSS_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024,resumed_from=args.resume_from,measurement_scope='This successful driver invocation only; failed presentation renders and the short render smoke check are not included. Earlier completed analysis and figures have separate resource records.'));print('POSTPROCESS_DONE',flush=True)
