import json
from pathlib import Path
import pytest
from PIL import Image, ImageStat
from fem3d.audit import sha256

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("name",["boundary_tags.png","tetrahedral_cutaway.png","mesh_quality_distribution.png","port_closeups.png","worst_elements.png"])
def test_review_image_is_nonempty_and_bound_to_mesh(name):
    base=ROOT/"reports/stage01"
    manifest=json.loads((base/"visualization_manifest.json").read_text())
    path=base/name
    with Image.open(path) as im:
        assert im.width>=1000 and im.height>=700
        assert max(ImageStat.Stat(im.convert("RGB")).stddev)>10
    assert path.stat().st_size>10000
    assert sha256(path)==manifest["images"][name]["sha256"]
    assert manifest["flow_fields_created"] is False
    assert manifest["cutaway"]["slice_polygons"]>0 and manifest["cutaway"]["cut_cells"]>0
