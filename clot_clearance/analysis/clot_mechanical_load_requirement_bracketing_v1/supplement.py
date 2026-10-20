"""Small evidence-based geometric refinement, then two normal mechanisms.

No random load search, no new shape/area/duration/material parameters.
"""
from common import *
from campaign_recovered import execute

def summaries(mode):
 out=[]
 for p in (O/'cases').glob('*/SUMMARY.json'):
  z=json.loads(p.read_text())
  if z['params']['mode']==mode:out.append(z)
 return out

def main():
 assert (O/'RECOVERED_DECADE_SCAN_FINISHED.json').exists()
 decisions=[]
 for iteration in range(2):
  cc=summaries('SURFACE_SHEAR');pos=[z['params']['peak_pa'] for z in cc if 'first_softening' in z['events']]
  if not pos:break
  hi=min(pos);neg=[z['params']['peak_pa'] for z in cc if z['status']=='COMPLETED' and 'first_softening' not in z['events'] and z['params']['peak_pa']<hi]
  if not neg:break
  lo=max(neg)
  # A completed no-event case above a positive contradicts a monotone onset bracket.
  if any(z['status']=='COMPLETED' and 'first_softening' not in z['events'] and z['params']['peak_pa']>hi for z in cc):break
  peak=float(np.sqrt(lo*hi));decisions.append(dict(mode='SURFACE_SHEAR',peak_pa=peak,reason='Geometric midpoint of observed monotone onset bracket',lower_pa=lo,upper_pa=hi));dump(O/'SUPPLEMENT_DECISIONS.json',decisions)
  z=execute('SURFACE_SHEAR',peak)
  if z['status']!='COMPLETED':break
 cc=summaries('SURFACE_SHEAR');completed=[z for z in cc if z['status']=='COMPLETED'];damage=[z for z in completed if 'first_softening' in z['events']]
 assert damage,'Manual scientific review required before adding normal-mode levels without a shear damage scale'
 # Preserve common amplitudes across mechanisms for an interpretable first contrast.
 levels=[100.,1000.]
 # Register the independent normal-mode requests serially before parallel launch.
 # execute() will find these existing entries, so there is no shared-book writer race.
 import concurrent.futures,time
 book=json.loads((O/'LOAD_CASES.json').read_text())
 for mode in ['LOCAL_IMPINGEMENT','IDEAL_LOCAL_PULL']:
  for peak in levels:
   tag=mode.lower()+'_'+format(peak,'.8g').replace('.','p')+'pa_r2'
   if not any(q['tag']==tag for q in book['cases']):book['cases'].append(dict(tag=tag,mode=mode,peak_pa=peak,request=O/'requests'/(tag+'.json'),registration='PLANNED; executed only if SUMMARY exists'))
   decisions.append(dict(mode=mode,peak_pa=peak,reason='Two common decade levels spanning observed shear damage; normal modes independent and each case uses one rank/thread'))
 dump(O/'LOAD_CASES.json',book);dump(O/'SUPPLEMENT_DECISIONS.json',decisions)
 def normal_mode(mode):
  begin=time.time();out=[]
  for peak in levels:
   z=execute(mode,peak);out.append(dict(peak_pa=peak,status=z['status']))
   if z['status']!='COMPLETED' or any(len(e['parent_atom_ids'])>=14 for e in z['secondary_fracture_events']):break
  return dict(mode=mode,start_unix_s=begin,end_unix_s=time.time(),cases=out)
 begin=time.time()
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:concurrency=list(pool.map(normal_mode,['LOCAL_IMPINGEMENT','IDEAL_LOCAL_PULL']))
 dump(O/'NORMAL_MODE_CONCURRENCY.json',dict(max_simultaneous_cases=2,each_case_ranks=1,each_case_threads=1,start_unix_s=begin,end_unix_s=time.time(),mode_workers=concurrency))
 dump(O/'SUPPLEMENT_FINISHED.json',dict(status='REVIEW_BEFORE_ANY_EXTRA_CASE',decisions=decisions))
if __name__=='__main__':main()
