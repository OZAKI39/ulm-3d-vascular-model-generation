from sv13j_support import *
from PIL import Image
def test_ten_images_and_hashes():
    d=load('visuals');assert len(d['figures'])==10
    for f in d['figures']:
        p=ROOT/f['path'];assert sha256(p)==f['sha256']
        with Image.open(p) as im:assert im.width>=1200 and im.height>=700;im.verify()
def test_unknown_scientific_panels_explicitly_not_run():
    d=load('visuals')
    panels=[f for f in d['figures'] if f['label']=='NOT RUN']
    assert len(panels)==5 and all(not f['has_measurements'] and f['numeric_bars']==0 for f in panels)

