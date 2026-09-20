from sv13l_support import *
def test_required_actual_figures():
    d=load('visuals');assert len(d['figures'])>=4
    for f in d['figures']:
        p=ROOT/f['path'];assert p.stat().st_size>10000 and sha256(p)==f['sha256']
