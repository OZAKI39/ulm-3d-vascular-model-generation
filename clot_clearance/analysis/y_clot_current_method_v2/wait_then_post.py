from common import *
import time,subprocess,sys,datetime

summary=O/'cases/biphasic/SUMMARY.json'
print('WAITING_FOR_THE_SINGLE_FORMAL_RUN',datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
while not summary.exists():time.sleep(10)
dump(O/'POSTPROCESS_SOURCE_FREEZE.json',{str(p.relative_to(R)):sha(p) for p in
    [A/'postprocess.py',A/'finish.py',R/'analysis/clot_fragment_surface_coupling_v3_red_visuals/render.py',R/'analysis/clot_fragment_surface_coupling_v3/render_videos.py']})
for name in ['postprocess.py','finish.py']:
    subprocess.run([sys.executable,'-u',str(A/name)],check=True)
print('POSTPROCESS_AND_REPORT_COMPLETE',datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)

