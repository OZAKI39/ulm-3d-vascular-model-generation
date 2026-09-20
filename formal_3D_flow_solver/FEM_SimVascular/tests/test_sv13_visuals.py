from sv13_support import artifact,ROOT
from sv_validation.provenance import sha256
from PIL import Image
def test_ten_review_figures_and_provenance():
    d=artifact('visuals');assert len(d['figures'])==10
    for f in d['figures']:
        p=ROOT/f['path'];assert sha256(p)==f['sha256']
        with Image.open(p) as im:assert im.width>=1000 and im.height>=600
    assert d['unavailable_measurements_not_plotted_as_zero']
