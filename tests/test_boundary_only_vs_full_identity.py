from __future__ import annotations

import json
from pathlib import Path


def test_boundary_only_metric_is_retained_but_not_final_referee() -> None:
    root = Path(__file__).resolve().parents[1]
    qc = (
        root
        / "outputs"
        / "cfd_flow"
        / "healthy_mouse_capillary_tau1_base_anchor003274_20260830"
        / "qc"
    )
    replay = json.loads(
        (qc / "tau1_full_timestep_replay_8step.json").read_text(encoding="utf-8")
    )
    final = json.loads((qc / "tau1_referee_v2_final.json").read_text(encoding="utf-8"))

    assert replay["R_boundary_only_max"] > replay["R_full_one_step_identity_max"]
    assert final["boundary_write_only_is_diagnostic_not_final_identity"] is True
    assert final["identity_scope"] == "complete source-proven Musubi timestep"
