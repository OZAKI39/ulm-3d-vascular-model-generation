from common import *
import subprocess,sys

def execute(mode,peak):
 tag=mode.lower()+'_'+format(peak,'.8g').replace('.','p')+'pa';case=O/'cases'/tag;request=O/'requests'/(tag+'.json');dump(request,dict(case=case,params=dict(mode=mode,peak_pa=peak)))
 book=json.loads((O/'LOAD_CASES.json').read_text())
 if not any(q['tag']==tag for q in book['cases']):
  book['cases'].append(dict(tag=tag,mode=mode,peak_pa=peak,request=request));dump(O/'LOAD_CASES.json',book)
 rc=subprocess.run([sys.executable,str(A/'native_case.py'),'--request',str(request)]).returncode
 if not (case/'SUMMARY.json').exists():raise RuntimeError(f'No summary for {tag}, rc={rc}')
 return json.loads((case/'SUMMARY.json').read_text())
if __name__=='__main__':
 for p in [1.,10.,100.,1000.]:
  z=execute('SURFACE_SHEAR',p)
  if z['status']!='COMPLETED' or z['classification']=='SECONDARY_FRAGMENTATION':break
 dump(O/'MAIN_DECADE_SCAN_FINISHED.json',dict(status='READY_FOR_EVIDENCE_BASED_REFINEMENT'))
