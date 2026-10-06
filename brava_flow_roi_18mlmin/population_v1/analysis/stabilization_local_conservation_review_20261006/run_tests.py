from common import *
import inspect,time,traceback
import test_review
results=[]
for name,fn in inspect.getmembers(test_review,inspect.isfunction):
    if not name.startswith('test_'):continue
    t=time.monotonic()
    try:fn();status='PASS';error=None
    except Exception:status='FAIL';error=traceback.format_exc()
    results.append(dict(test=name,status=status,elapsed_s=time.monotonic()-t,error=error));print(name,status,error or '',flush=True)
put('TEST_RESULTS.json',dict(status='PASS' if all(r['status']=='PASS' for r in results) else 'FAIL',count=len(results),tests=results,scope='Engineering consistency only; no physical certification.',patch_tests='NOT_APPLICABLE: solver not patched'))
if any(r['status']!='PASS' for r in results):raise SystemExit(1)
