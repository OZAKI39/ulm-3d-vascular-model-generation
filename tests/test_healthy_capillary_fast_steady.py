from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from utils.cfd_flow.healthy_capillary_calibration import OUTLET_GAUGE_PRESSURE_PA
from utils.cfd_flow.healthy_capillary_fast_steady import (
    BOUNDARY_HARD,
    LONG_RANGE,
    MASS_HARD,
    PRESSURE_HARD,
    SHORT_RANGE,
    SOURCE_ITERATION,
    TARGET_MASS_FLOW_KG_S,
    VELOCITY_HARD,
    FlushPolicy,
    RestartEvidence,
    all_final_gates_pass,
    boundary_residual,
    characteristic_pressure_drop,
    choose_mpi_ranks,
    classify_stage,
    create_stop_file,
    inlet_residual,
    mass_residual,
    next_benchmark_rank,
    pressure_residual,
    restart_compatibility,
    select_window,
    significant_backflow,
    velocity_l2_residual,
)


def _evidence(iteration: int) -> RestartEvidence:
    path = Path(f"checkpoint_{iteration}")
    return RestartEvidence(iteration, path.with_suffix(".lua"), path.with_suffix(".lsb"), "sha", "test")


def test_mpi_candidate_selection_and_expansion() -> None:
    assert next_benchmark_rank([], 12) == 4
    rows = [{"ranks": 4, "iterations_per_s": 20.0}, {"ranks": 6, "iterations_per_s": 30.0}]
    assert next_benchmark_rank(rows, 12) == 8
    rows.append({"ranks": 8, "iterations_per_s": 32.0})
    assert next_benchmark_rank(rows, 12) == 12
    rows[-1]["iterations_per_s"] = 29.0
    assert next_benchmark_rank(rows, 12) is None


def test_rank_tie_rule_under_three_percent_chooses_fewer() -> None:
    result = choose_mpi_ranks([
        {"ranks": 4, "iterations_per_s": 100.0},
        {"ranks": 6, "iterations_per_s": 102.0},
        {"ranks": 8, "iterations_per_s": 101.0},
    ])
    assert result["selected"]["ranks"] == 4


def test_restart_compatibility_for_frozen_resume() -> None:
    root = Path(__file__).resolve().parents[1]
    restart = root / "outputs" / "cfd_flow" / "healthy_mouse_capillary_calibration_anchor003274_20260829_180310" / "restart"
    result = restart_compatibility(restart / "a3274_lastHeader.lua", restart / "a3274_7.670E-03.lsb")
    assert result["status"] == "PASS"
    assert result["iteration"] == SOURCE_ITERATION


def test_short_and_long_window_selectors() -> None:
    candidates = [_evidence(299_326), _evidence(309_326), _evidence(314_166)]
    short = select_window(319_166, candidates, target=10_000, allowed=SHORT_RANGE)
    long = select_window(319_166, candidates, target=20_000, allowed=LONG_RANGE)
    assert short is not None and short.iteration == 309_326
    assert long is not None and long.iteration == 299_326
    assert select_window(322_000, candidates, target=10_000, allowed=SHORT_RANGE) is None


def test_mass_residual_physical_conversion() -> None:
    from utils.cfd_flow.healthy_capillary_fast_steady import EXPECTED_DT_S, EXPECTED_DX_M, REFERENCE_DENSITY_KG_M3
    delta_iter = 10_000
    desired_accumulation = 0.005 * TARGET_MASS_FLOW_KG_S
    delta_total = desired_accumulation * delta_iter * EXPECTED_DT_S / (REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**3)
    assert math.isclose(mass_residual(100.0, 100.0 + delta_total, delta_iter), 0.005, rel_tol=1e-10)


def test_velocity_l2_residual_uses_full_field() -> None:
    previous = np.zeros((2, 3))
    current = np.asarray(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)))
    assert math.isclose(velocity_l2_residual(previous, current), 1.0)
    assert velocity_l2_residual(current * 0.995, current) <= VELOCITY_HARD


def test_pressure_characteristic_drop_and_residual() -> None:
    inlet = 301.0
    expected = float(np.median([abs(inlet - value) for value in OUTLET_GAUGE_PRESSURE_PA.values()]))
    characteristic = characteristic_pressure_drop(inlet)
    assert characteristic == expected
    assert pressure_residual(0.004 * characteristic, characteristic) <= PRESSURE_HARD


def test_inlet_boundary_and_backflow_rules() -> None:
    inlet = TARGET_MASS_FLOW_KG_S
    assert inlet_residual(inlet * 1.005) < 0.01
    outlets = [0.2 * inlet, 0.5 * inlet, 0.3 * inlet]
    assert boundary_residual(inlet, outlets) < BOUNDARY_HARD
    assert not significant_backflow([0.5 * inlet, -0.04 * inlet, 0.54 * inlet], inlet)
    assert significant_backflow([0.5 * inlet, -0.051 * inlet, 0.551 * inlet], inlet)


def test_stage_and_all_gate_classification() -> None:
    assert classify_stage(0.051) == "STAGE_1_FAR"
    assert classify_stage(0.02) == "STAGE_2_NEAR"
    assert classify_stage(0.01) == "STAGE_3_FINAL"
    record = {
        "R_mass_short": MASS_HARD, "R_mass_long": MASS_HARD,
        "R_boundary": BOUNDARY_HARD, "R_velocity": VELOCITY_HARD,
        "R_pressure": PRESSURE_HARD, "R_inlet": 0.01,
        "pressure_pdf_crosscheck": {"relative_discrepancy": 0.001, "direction_same": True},
        "significant_backflow": False, "all_finite": True, "minimum_pdf": 0.01,
        "maximum_lattice_speed": 0.001, "inlet_globbc": 287,
    }
    assert all_final_gates_pass(record)
    record["R_boundary"] = BOUNDARY_HARD + 1e-8
    assert not all_final_gates_pass(record)


def test_stop_file_and_gate_pass_checkpoint_selection(tmp_path: Path) -> None:
    path = create_stop_file(tmp_path)
    assert path == tmp_path / "stop" and path.is_file()
    candidates = [_evidence(309_326), _evidence(314_166)]
    selected = select_window(319_166, candidates, target=10_000, allowed=SHORT_RANGE)
    assert selected is not None and selected.iteration == 309_326


class _FakeStream:
    def __init__(self) -> None:
        self.flushes = 0

    def flush(self) -> None:
        self.flushes += 1


def test_writer_buffering_count_or_time() -> None:
    now = [0.0]
    stream = _FakeStream()
    policy = FlushPolicy(stream, record_limit=3, interval_s=1.0, clock=lambda: now[0])
    policy.note()
    policy.note()
    assert stream.flushes == 0
    policy.note()
    assert stream.flushes == 1
    now[0] = 1.1
    policy.note()
    assert stream.flushes == 2
