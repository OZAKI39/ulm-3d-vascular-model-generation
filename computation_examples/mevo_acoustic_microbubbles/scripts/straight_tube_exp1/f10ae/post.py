"""Bounded scientific review and artifacts after full frozen-action campaign."""
import time
from .common import *
def main():
    begin=time.perf_counter()
    while not (BASE/'WORKFLOW_COMPLETE.json').exists():
        if time.perf_counter()-begin>1800:raise RuntimeError('Scientific workflow did not finish in allotted wait')
        time.sleep(3)
    from . import low_admit,audit,interpret,figures,videos,report,finalize,diagnostics
    for name,module,marker in [('low_admit',low_admit,'LOW_PRESSURE_ADMISSION.json'),('audit',audit,'AUDIT.json'),('diagnostics',diagnostics,'ROBUSTNESS_FIELD_FORCE.json'),('interpret',interpret,'INTERPRETATION.json'),('figures',figures,'figures/MANIFEST.json'),('videos',videos,'videos/MANIFEST.json'),('report',report,'FINAL_CHINESE_REPORT.md')]:
        if (BASE/marker).exists():continue
        write(BASE/'POST_PROGRESS.json',dict(phase=name,identity=identity(),elapsed_s=time.perf_counter()-begin));module.main()
    finalize.main()
    print('POST_COMPLETE',time.perf_counter()-begin,flush=True)
if __name__=='__main__':main()
