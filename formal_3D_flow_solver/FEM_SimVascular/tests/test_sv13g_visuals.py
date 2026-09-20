from sv13g_support import *
from sv_validation.provenance import sha256
def test_eight_figures_and_honest_not_run_panels():
    d=load('visuals');assert len(d['figures'])==8
    for f in d['figures']:
        p=ROOT/f['path'];assert p.stat().st_size>1000 and sha256(p)==f['sha256']
        if not f['has_measurements']:assert f['label']=='NOT RUN' and f['numeric_bars']==0
