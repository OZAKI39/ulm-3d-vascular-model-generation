import json
import pytest
from PIL import Image,ImageStat
from stage02_helpers import ROOT
from fem3d.audit import sha256


@pytest.mark.parametrize("name",["pipe_geometry_and_bc.png","velocity_profile.png","pressure_axial.png","lambda_vs_Q.png","mesh_convergence.png","flux_balance.png","velocity_slice_3d.png"])
def test_review_figure_is_readable_and_nonempty(name):
    base=ROOT/"reports/stage02"
    manifest=json.loads((base/"visualization_manifest.json").read_text())
    with Image.open(base/name) as im:
        assert im.width>=1000 and im.height>=650
        assert max(ImageStat.Stat(im.convert("RGB")).stddev)>10
    assert sha256(base/name)==manifest["images"][name]["sha256"]
