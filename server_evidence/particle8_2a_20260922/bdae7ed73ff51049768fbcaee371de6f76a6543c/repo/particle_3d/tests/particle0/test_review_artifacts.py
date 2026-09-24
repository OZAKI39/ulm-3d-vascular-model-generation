"""Keep every requested review artifact and its plotted data reproducible."""
from pathlib import Path
import csv
import json
import numpy as np
from PIL import Image
import pytest

REPORT = Path(__file__).resolve().parents[2] / "reports/particle0"
PAIRS = [
    ("00_frozen_input_overview.png", "00_frozen_input_summary.json"),
    ("01_affine_field_validation.png", "01_affine_field_validation.csv"),
    ("02_node_back_sampling_error.png", "02_node_back_sampling_error.csv"),
    ("03_shared_face_continuity.png", "03_shared_face_continuity.csv"),
    ("04_inside_outside_classification.png", "04_inside_outside_classification.csv"),
    ("05_real_flow_velocity_vectors.png", "05_real_flow_samples.csv"),
    ("06_real_flow_scalar_diagnostics.png", "05_real_flow_samples.csv"),
    ("07_sampling_runtime_smoke.png", "07_sampling_runtime_smoke.csv"),
]


@pytest.mark.parametrize("figure,data", PAIRS)
def test_permanent_screen_readable_figure_and_nonempty_data(figure, data):
    with Image.open(REPORT / "figures" / figure) as im:
        assert im.format == "PNG" and im.width >= 1500 and im.height >= 1000
        assert np.std(np.asarray(im)) > 10  # reject blank images
    path = REPORT / "data" / data
    if path.suffix == ".json":
        assert json.loads(path.read_text())
    else:
        assert b"\r" not in path.read_bytes()  # keep repository CSV line endings LF
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) >= 3


def test_all_seven_steps_have_chinese_note_and_numerical_metrics():
    for stage in range(7):
        notes = list(REPORT.glob(f"{stage:02d}_*.md"))
        assert len(notes) == 1
        assert "检查说明" in notes[0].read_text()
        metrics = json.loads((REPORT / "data" / f"metrics_p0{stage}.json").read_text())
        assert metrics.get("passed", metrics.get("status") == "PASS")
