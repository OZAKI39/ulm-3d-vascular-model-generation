"""Numerical Jacobian equivalence, frozen acceptance, and bounded fitting timing."""
import sys,time,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from fit_adapter import fit_affine as fast
from scipy.optimize._numdiff import approx_derivative
import cv2
rng=np.random.default_rng(200073);rows=[]
for i in range(80):
 n=int(rng.integers(8,90));xy=rng.uniform(24,216,(n,2));X=np.c_[np.ones(n),(xy-120)/120];coef=rng.normal(0,.3,(3,2));d=X@coef+rng.normal(0,[.02,.4,2.][i%3],(n,2));start=rng.normal(size=6)
 analytic=np.kron(X,np.eye(2));numeric=approx_derivative(lambda a:(X@a.reshape(3,2)-d).ravel(),start)
 assert np.max(np.abs(analytic-numeric))<1e-8
 t=time.perf_counter();a=fit_affine(xy,d,.0999,(240,240),'frozen');old_ms=(time.perf_counter()-t)*1000
 t=time.perf_counter();b=fast(xy,d,.0999,(240,240),'exact_jac');new_ms=(time.perf_counter()-t)*1000
 diff=np.max(np.abs(np.asarray(a.get('raw_estimate',[0]*6))-np.asarray(b.get('raw_estimate',[0]*6))))
 rows.append(dict(index=i,old_ms=old_ms,new_ms=new_ms,max_six_difference=float(diff),same_valid=a['valid']==b['valid'],same_reasons=a['reasons']==b['reasons'],old_residual=a['residual'],new_residual=b['residual']))
 assert diff<1e-4 and a['valid']==b['valid'] and a['reasons']==b['reasons'],rows[-1]
dump(OUT/'results/fit_adapter_check.json',dict(identity=identity(),rows=rows,old_ms=stats([x['old_ms'] for x in rows]),new_ms=stats([x['new_ms'] for x in rows]),max_six_difference=max(x['max_six_difference'] for x in rows),all_validity_equal=True,scope='Identical F10U bytecode/gates/logm with SciPy exact derivative. No FFT or scientific model rewrite.'))
print('fit equivalence passed; old/new p95',stats([x['old_ms'] for x in rows]),stats([x['new_ms'] for x in rows]))
