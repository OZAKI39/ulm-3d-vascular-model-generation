"""Start the already-defined postprocess once after the existing formal service exits."""
from common import *
import subprocess,time
FORMAL='y_clot_irregular_v4_20261020T041100Z'
POST='y_clot_v4_postprocess_20261020T041100Z'
while True:
 status=subprocess.run(['supervisorctl','status',FORMAL],capture_output=True,text=True).stdout.strip()
 if ' RUNNING ' in status or ' STARTING ' in status:
  time.sleep(45);continue
 result=json.loads((O/'NATIVE_V4_RESULT.json').read_text())
 assert result['status'] in ['COMPLETED','NUMERICAL_FAILURE'],('Formal service stopped without terminal result',status,result['status'])
 assert not (O/'POSTPROCESS_AUTO_START.json').exists(),'Postprocess start already requested'
 dump(O/'POSTPROCESS_AUTO_START.json',dict(formal_service_status=status,result_status=result['status'],one_request=True))
 subprocess.run(['supervisorctl','start',POST],check=True)
 print('POSTPROCESS_STARTED_AFTER_TERMINAL_FORMAL_RESULT',flush=True)
 break
