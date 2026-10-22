"""Serialized remaining physical experiments after independent validation."""
import time,traceback
from .common import *
def main():
    begin=time.perf_counter();write(BASE/'REMAINING_WORKFLOW_START.json',dict(identity=identity()))
    deadline=begin+1500
    while not (BASE/'GENERALIZATION_RESULTS.json').exists():
        if time.perf_counter()>deadline:raise RuntimeError('Independent validation did not complete in its registered budget')
        time.sleep(3)
    from .longtime import main as longtime
    from .sobol import main as sobol
    from .audit import main as audit
    for name,fn in [('longtime',longtime),('sobol',sobol),('audit',audit)]:
        try:fn()
        except Exception as ex:
            write(BASE/(name.upper()+'_WORKFLOW_FAILURE.json'),dict(reason=repr(ex),traceback=traceback.format_exc(),identity=identity()));raise
    write(BASE/'PHYSICS_WORKFLOW_COMPLETE.json',dict(wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
