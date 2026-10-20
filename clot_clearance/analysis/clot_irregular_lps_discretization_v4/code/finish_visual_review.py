"""Finish reviewed visual changes only. Never start or repeat dynamics."""
from common import *
import subprocess,sys,datetime

def main():
 result=json.loads((O/'NATIVE_V4_RESULT.json').read_text());assert result['status']=='COMPLETED'
 physical_sha=sha(O/'NATIVE_V4_RESULT.json');record=dict(status='VISUAL_REVIEW_REFRESH',native_result_sha256=physical_sha,steps=[])
 dump(O/'POSTPROCESS_STATUS.json',record)
 try:
  for args in [['analyze_result.py','--spectrum-only'],['render_videos.py'],['render_target_closeup.py'],['write_report.py'],['validate_delivery.py']]:
   print('VISUAL_REVIEW',args,flush=True);subprocess.run([sys.executable,str(A/args[0]),*args[1:]],check=True,cwd=R);record['steps'].append(args)
  assert sha(O/'NATIVE_V4_RESULT.json')==physical_sha
  record.update(status='COMPLETE',accepted_time_s=result['accepted_time_s'],finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),no_dynamics_rerun=True)
 except BaseException as exc:
  record.update(status='FAILED_VISUAL_REVIEW',exception=repr(exc));raise
 finally:dump(O/'POSTPROCESS_STATUS.json',record)
 print('VISUAL_REVIEW_COMPLETE',flush=True)

if __name__=='__main__':main()
