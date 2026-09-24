from PIL import Image,ImageStat
from fem3d.audit import sha256
from stage015_helpers import REPORT,read


def test_review_images_have_source_manifest_and_clear_failure_labels():
    manifest=read(REPORT/'visualization_manifest.json')
    assert len(manifest['images'])==8
    decision=read(REPORT/'candidate_selection.json')
    for name,item in manifest['images'].items():
        assert sha256(REPORT/name)==item['sha256']
        with Image.open(REPORT/name) as im:
            assert im.width>=1000 and im.height>=650
            assert max(ImageStat.Stat(im.convert('RGB')).stddev)>10
        if decision['selected_candidate'] is None: assert item['selected_candidate'] is None
