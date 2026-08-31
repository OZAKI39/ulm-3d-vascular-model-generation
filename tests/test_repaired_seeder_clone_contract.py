from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.tau1_grid_convergence import (
    BASE_MESH_RUN,
    GRID_SPECS,
    render_repaired_seeder_config,
    seeder_physical_spatial_signature,
)


ROOT = Path(__file__).resolve().parents[1]


def test_coarse_and_fine_change_only_lattice_geometry() -> None:
    path = ROOT / "outputs" / "cfd_flow" / BASE_MESH_RUN / "seeder" / "seeder.lua"
    base = path.read_text(encoding="utf-8")
    signature = seeder_physical_spatial_signature(base)
    for label in ("coarse", "fine"):
        rendered = render_repaired_seeder_config(base, GRID_SPECS[label])
        assert seeder_physical_spatial_signature(rendered) == signature
        assert "calc_dist = true" in rendered
