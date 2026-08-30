import pytest

from utils.cfd_flow.musubi_wall_force_diagnostics import wall_force_decision


@pytest.mark.parametrize(
    ("arguments", "expected"),
    [
        ({"force_pass": False, "wall_pass": None, "tau1_pass": None}, "STOP_FORCE_CONTRACT"),
        ({"force_pass": True, "wall_pass": True, "tau1_pass": None}, "RUN_TAU1_N27"),
        (
            {"force_pass": True, "wall_pass": True, "tau1_pass": True},
            "HIGH_TAU_BGK_WALL_COUPLING_CONFIRMED",
        ),
        (
            {
                "force_pass": True,
                "wall_pass": True,
                "tau1_pass": False,
                "official_pass": False,
            },
            "CFD_FLOW_UPSTREAM_PIPE_FORCE_REFERENCE_FAILED",
        ),
    ],
)
def test_decision_tree_has_no_unrequested_sweep(arguments: dict[str, object], expected: str) -> None:
    assert wall_force_decision(**arguments) == expected
