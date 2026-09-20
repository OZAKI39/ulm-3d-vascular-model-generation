"""Bounded sequential matrix: four early smokes, healthy late windows, only top-two early windows."""
import json,subprocess,sys
from remote import ROOT,upload
S=ROOT/'scripts/sv13q';R=ROOT/'reports/sv1_3q'
def read(n):return json.loads((R/(n+'.json')).read_text())
def write(n,d):(R/(n+'.json')).write_text(json.dumps(d,indent=2)+'\n')
def run(c,mode):
 name=c+'_'+mode.upper()
 assert not (R/(name+'_acceptance.json')).exists(),'No implicit CFD repeats'
 p=subprocess.run([sys.executable,'-B',S/'run_case.py',c,mode])
 assert (R/(name+'_acceptance.json')).exists(),'Engineering failure: inspect and repair before continuing'
 d=read(name+'_acceptance');print(c,mode,d['status'],d['wall_time_s'],flush=True)
 return d
assert read('remote/svmp_reuse_build')['status']=='PASS'
for n in ('flow_parser.py',):upload(S/n,n)
for c in ('R2','R3','R5','RA'):run(c,'smoke')
late=[]
for c in ('R2','R3','R5','RA'):
 if read(c+'_SMOKE_acceptance')['status']=='PASS':
  d=run(c,'window')
  if d['status']=='PASS':late.append(d)
late.sort(key=lambda d:d['wall_time_s']);top=[d['candidate'] for d in late[:2]]
write('late_ranking',dict(status='PASS',ranking=[dict(candidate=d['candidate'],wall_time_s=d['wall_time_s']) for d in late],top_two=top,selection='Only healthy complete ten-step late windows'))
for c in top:run(c,'early')
subprocess.run([sys.executable,'-B',S/'select_winner.py'],check=True)
print('All authorized windows completed. Full run requires separate >=10% gate.',flush=True)
