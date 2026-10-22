"""Final read-only diagnostic/report revision, after the first artifact pass."""
import time,shutil
from .common import *
def main():
    start=time.perf_counter()
    while not (BASE/'INPUT_OUTPUT_SHA256.json').exists():
        if time.perf_counter()-start>1800:raise RuntimeError('First artifact pass incomplete')
        time.sleep(3)
    from . import diagnostics,report,figures,finalize
    if not (BASE/'ROBUSTNESS_FIELD_FORCE.json').exists():diagnostics.main()
    prior=BASE/'evidence/report_before_final_review';prior.mkdir(exist_ok=False)
    for p in BASE.glob('*.md'):shutil.copy2(p,prior/p.name)
    # Readability correction for the non-extending OFF reference; no data altered.
    for p in (BASE/'figures').iterdir():
        if p.is_file():shutil.copy2(p,prior/p.name)
    report.main();figures.main();finalize.main()
    print('REVIEW_COMPLETE',time.perf_counter()-start,flush=True)
if __name__=='__main__':main()
