"""Stage Q tests inspect immutable artifacts; never invoke CFD or SSH."""
import hashlib,json,re
from pathlib import Path
from sv_validation.sv13q import *
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'reports/sv1_3q';C=ROOT/'configs/sv1_3q'
def read(name):return json.loads((R/(name+'.json')).read_text())
def policy():return json.loads((C/'policy.json').read_text())
def case(candidate,mode='WINDOW'):return read(candidate+'_'+mode+'_acceptance')
def cases():return [json.loads(p.read_text()) for p in sorted(R.glob('*_acceptance.json'))]
def raw(d):return (ROOT/'logs/sv1_3q/remote'/(d['name']+'.log')).read_text()
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check_candidate(c):
 smoke=case(c,'SMOKE');assert smoke['steps']<=3 and smoke['start_step']==0
 if smoke['status']=='FAIL':assert smoke['errors'];assert not (R/(c+'_WINDOW_acceptance.json')).exists();return
 late=case(c);assert late['start_step']==60
 if late['status']=='PASS':
  assert late['steps']==10 and late['stop_step']==70
  assert late['reuse']['status']=='PASS'
 else:assert late['errors']
