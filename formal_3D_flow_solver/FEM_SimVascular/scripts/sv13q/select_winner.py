import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3q'
def read(n):return json.loads((R/(n+'.json')).read_text())
base=read('baseline_r1');total=base['late']['wall_time_s']+base['early']['wall_time_s'];eligible=[];rejected=[]
for c in read('late_ranking')['top_two']:
 early=read(c+'_EARLY_acceptance');late=read(c+'_WINDOW_acceptance')
 reason='EARLY_UNHEALTHY' if early['status']!='PASS' else 'EARLY_OBSERVED_MORE_THAN_5_PERCENT_SLOWER_THAN_R1' if early['wall_time_s']>1.05*base['early']['wall_time_s'] else None
 if reason:rejected.append(dict(candidate=c,reason=reason));continue
 wall=early['wall_time_s']+late['wall_time_s'];eligible.append(dict(candidate=c,observed_total_wall_s=wall,late_wall_s=late['wall_time_s'],early_wall_s=early['wall_time_s'],improvement=1-wall/total,PETSC_OPTIONS=early['PETSC_OPTIONS'],build_report='svmp_reuse_build'))
eligible.sort(key=lambda d:d['observed_total_wall_s']);best=eligible[0] if eligible else {}
d=dict(status='SELECTED' if best else 'NO_ELIGIBLE_WINNER',baseline_total_wall_s=total,eligible=eligible,rejected=rejected,full_required=bool(best and best['improvement']>=.10),maximum_full_runs=1,observational=True,early_baseline_limitation=base['early']['limitation'],**best)
(R/'winner.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d),flush=True)
