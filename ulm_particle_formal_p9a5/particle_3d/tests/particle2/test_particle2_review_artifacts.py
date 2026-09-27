import json
from pathlib import Path
import numpy as np
import pytest
from PIL import Image
from particle_3d.rbc_distribution import digest


@pytest.mark.parametrize("index",range(11))
def test_each_figure_has_readable_png_bound_source_data_and_chinese_note(index,p2_data):
    source=json.loads((p2_data/f"{index:02d}_figure_sources.json").read_text())
    with Image.open(p2_data.parent/"figures"/source["figure"]) as image:
        assert image.format=="PNG" and image.width>=1500 and image.height>=900
        assert np.std(np.asarray(image))>10
    for path,expected in source["source_sha256"].items():
        assert digest(p2_data/path)==expected
    note=p2_data.parent/("00_distribution_literature_report.md" if index==0 else f"{index:02d}_step_notes.md")
    assert note.is_file() and any("\u4e00"<=c<="\u9fff" for c in note.read_text())
    assert source["manual_visual_review"]=="PENDING_USER_REVIEW"


def test_population_metadata_binds_csv_and_contract(p2_data,population,p2_repo):
    meta=json.loads((p2_data/"C57BL6_RBC_GEOMETRY_VALIDATION_100000.metadata.json").read_text())
    assert meta["N"]==100000 and meta["seed"]==2026092002
    assert meta["csv_sha256"]==digest(p2_data/meta["csv_file"])
    assert meta["contract_sha256"]==digest(p2_repo/"particle_3d/contracts/C57BL6_RBC_GEOMETRY_DISTRIBUTION_V0.json")
    assert meta["sample_structured_array_sha256"]==population.metadata["sample_structured_array_sha256"]
    assert not meta["formal_hematocrit_population"] and not meta["production_particle_population"]
