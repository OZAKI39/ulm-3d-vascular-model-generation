from analysis.dual_track_v1.common import *
import sys,time,traceback,resource
from analysis.dual_track_v1 import mesh_b,cfd_b
level=sys.argv[1];start=time.monotonic();stage='mesh';out=B/level
try:
 mesh_b.main(level);stage='CFD';save(out/'FORMAL_ATTEMPT.json',dict(formal_CFD_attempt=1,level=level,started_unix=time.time()))
 cfd_b.main(level);stage='COMPLETE'
except BaseException as e:
 save(B/(level+'_STOP.json'),dict(stage=stage,error=repr(e),traceback=traceback.format_exc(),wall_s=time.monotonic()-start));raise
finally:
 save(B/(level+'_EXECUTION.json'),dict(last_stage=stage,wall_s=time.monotonic()-start,peak_self_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,child_CPU_s=resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime))
