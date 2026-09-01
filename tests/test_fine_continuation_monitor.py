from __future__ import annotations

from utils.cfd_flow.fine_continuation_monitor import (
    CHECKPOINTS,
    EXPECTED_RANKS,
    TARGET_ITERATION,
    _classify,
    _next_checkpoint,
    restart_pair_complete,
)


def _probe(**updates: object) -> dict[str, object]:
    probe: dict[str, object] = {
        "fatal_matches": [],
        "ranks": [{"pid": value} for value in range(EXPECTED_RANKS)],
        "launcher_alive": True,
        "command_match": True,
        "runtime_status": "RUNNING",
        "successful_run": False,
        "latest_restart": {"iteration": 5000},
    }
    probe.update(updates)
    return probe


def test_healthy_running_requires_full_process_identity() -> None:
    assert _classify(_probe()) == "HEALTHY_RUNNING"
    assert _classify(_probe(command_match=False)) == "PROCESS_IDENTITY_MISMATCH"


def test_runtime_scientific_pattern_has_priority() -> None:
    assert _classify(_probe(fatal_matches=["nan"])) == "SCIENTIFIC_FAILURE"


def test_clean_completion_requires_target_restart_and_success_marker() -> None:
    stopped = {
        "fatal_matches": [],
        "ranks": [],
        "launcher_alive": False,
        "command_match": False,
        "runtime_status": "PASS",
        "successful_run": True,
        "latest_restart": {"iteration": TARGET_ITERATION},
    }
    assert _classify(stopped) == "CLEAN_COMPLETION"
    stopped["latest_restart"] = {"iteration": CHECKPOINTS[-2]}
    assert _classify(stopped) == "RECOVERABLE_OPERATIONAL_ERROR"


def test_next_checkpoint_uses_only_complete_restart_progress() -> None:
    assert _next_checkpoint(5000) == CHECKPOINTS[0]
    assert _next_checkpoint(CHECKPOINTS[0]) == CHECKPOINTS[1]
    assert _next_checkpoint(TARGET_ITERATION) is None


def test_restart_completeness_requires_header_iteration_and_exact_payload() -> None:
    values = {
        "header_text": "time_point = { sim = 1.0, iter = 202379 }",
        "header_size": 899,
        "payload_size": 400_949 * 19 * 8,
        "iteration": 202_379,
    }
    assert restart_pair_complete(**values)
    assert not restart_pair_complete(**(values | {"payload_size": values["payload_size"] - 8}))
    assert not restart_pair_complete(**(values | {"iteration": 202_380}))
