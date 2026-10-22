import time
from .common import *
def main():
    start=time.perf_counter()
    while not (BASE/'ANALYSIS_COMPLETE.json').exists():
        if time.perf_counter()-start>1700:raise RuntimeError('Bounded final analysis wait expired')
        time.sleep(3)
    from .repair_failed_view import main as repair
    from .rotation_checks import main as check
    repair();check();write(BASE/'POST_ANALYSIS_COMPLETE.json',dict(identity=identity()))
if __name__=='__main__':main()
