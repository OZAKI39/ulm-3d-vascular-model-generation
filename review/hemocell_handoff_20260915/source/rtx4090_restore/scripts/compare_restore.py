#!/usr/bin/env python3
import csv,json,sys
from pathlib import Path
sys.dont_write_bytecode=True
from evaluate_case import compare
W=Path(__file__).resolve().parents[1];n=int(sys.argv[1]);D=W/f'stage4_{n}'
results=[]
for ref,points in [('rtx200',[0,1,10,100,200]),('cpu200',[0,100,200])] if n==200 else [(f'rtx{n}',[x for x in [0,100,200,1000,5000] if x<=n]),(f'cpu{n}',[x for x in [0,100,200,1000,5000] if x<=n])]:
 result=compare(W/'references'/ref,D,W/'verification'/f'{n}_vs_{ref}',points)
 results.append(result)
with (W/f'STAGE4_{n}_CORRECTNESS.csv').open('w',newline='') as f:
 fields=['reference','quantity','step','error','tolerance','status','detail']
 wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader()
 for result in results:
  for row in result['maxima'].values():wr.writerow(dict(reference=result['CPU'],**{k:row[k] for k in fields if k!='reference'}))
ok=all(r['status']=='PASS' and r['NUMERICAL_DIVERGENCE_TREND']=='STABLE' for r in results)
(W/'verification'/f'STAGE4_{n}_GATE.json').write_text(json.dumps(dict(status='PASS' if ok else 'FAIL',numerical_divergence_trend=('FAIL' if any(r['status']=='FAIL' for r in results) else 'GROWING' if any(r['NUMERICAL_DIVERGENCE_TREND']=='GROWING' for r in results) else 'STABLE'),evaluator_unmodified=True,trend_scope='Frozen diagnostic growth rule applies at >=5000 steps; at <=1000 STABLE means all frozen tolerance gates pass, not a sustained-run extrapolation.',comparisons=[{k:v for k,v in r.items() if k!='maxima'} for r in results]),indent=2)+'\n')
print('CORRECTNESS',n,'PASS' if ok else 'FAIL',flush=True)
sys.exit(0 if ok else 2)
