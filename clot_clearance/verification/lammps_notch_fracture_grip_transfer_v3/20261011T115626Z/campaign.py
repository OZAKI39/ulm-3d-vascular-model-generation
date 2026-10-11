from pathlib import Path
import subprocess,time,json,datetime,resource,os
base=Path(__file__).resolve().parent
cmd=['python','-u','analysis/notch_grip_transfer_v3/run.py',str(base/'M5C')]
start=time.time()
r=subprocess.run(cmd,cwd='/workspace/y_clot/project')
u=resource.getrusage(resource.RUSAGE_CHILDREN)
(base/'RUN_RESOURCE.json').write_text(json.dumps(dict(command=cmd,start_utc=datetime.datetime.fromtimestamp(start,datetime.timezone.utc).isoformat(),end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),wall_s=time.time()-start,exit_code=r.returncode,child_user_s=u.ru_utime,child_system_s=u.ru_stime,max_rss_KiB=u.ru_maxrss,mpi_ranks=1,threads=1),indent=2))
raise SystemExit(r.returncode)
