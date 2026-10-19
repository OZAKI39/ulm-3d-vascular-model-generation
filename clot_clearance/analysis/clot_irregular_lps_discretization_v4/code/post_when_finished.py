"""Durable postprocessing of the one trajectory; never launches dynamics."""
from common import *
import subprocess,time,sys,traceback
def main():
 service='clot_v4_formal_20261019T221000Z'
 while True:
  q=subprocess.run(['supervisorctl','status',service],capture_output=True,text=True).stdout.split()
  if len(q)>=2 and q[1] in ['EXITED','FATAL']:break
  time.sleep(30)
 result=json.loads((O/'NATIVE_V4_RESULT.json').read_text());assert result['status'] in ['COMPLETED','NUMERICAL_FAILURE'],'No terminal native result; retain evidence and do not invent an endpoint'
 for script in ['terminal_audit.py','analyze_result.py','render_videos.py','write_report.py','validate_delivery.py']:
  print('POSTPROCESS',script,flush=True);subprocess.run([sys.executable,'-u',str(A/script)],check=True,cwd=A)
 dump(O/'POSTPROCESS_STATUS.json',dict(status='COMPLETE',native_status=result['status'],accepted_time_s=result['accepted_time_s']));print('ALL_POSTPROCESS_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except BaseException:
  dump(O/'POSTPROCESS_STATUS.json',dict(status='FAILED',traceback=traceback.format_exc()));raise
