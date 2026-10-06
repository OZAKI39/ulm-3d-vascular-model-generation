"""Development V3: accelerated residual-checked semi-implicit prefixes."""
from scripts.stage5a_r2.run import run
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write

if __name__=='__main__':
    rows=[]
    for case in ['B2','B3']:
        prior=None
        for end in [1e-6,1e-5,1e-4]:
            p,s=run(case,end,method='S1',label='prefix_v3',budget=600,resume=prior)
            rows.append(s);write(BASE/'diagnosis/PREFIX_V3_SUMMARY.json',rows)
            if s['status']!='COMPLETE':break
            prior=p/'CHECKPOINT.npz'
