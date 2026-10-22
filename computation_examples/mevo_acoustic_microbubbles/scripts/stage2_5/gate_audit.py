"""Materialize evidence-backed gate decisions; missing inputs close downstream gates."""

import hashlib
import json
from pathlib import Path

from mevo_mb.validation.gates import ValidationGateManager, ValidationGateError
from common import atomic_json, provenance


def evidence(path):
    p = Path(path)
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def main():
    root = Path("results/stage2_5")
    pair = root / "pair_asymptotics/STATE.json"
    zhang = root / "zhang_forensic/STATE.json"
    table = root / "response_table/STATE.json"
    trans = root / "translation_map/STATE.json"
    z = json.loads(zhang.read_text())
    t = json.loads(table.read_text()) if table.exists() else {}
    gates = {
        "PAIR_REFERENCE_PASS": dict(
            status="UNRESOLVED",
            scopes=[],
            evidence=[evidence(zhang)],
            reason="Fig3 nonlinear discrepancy persists; no demonstrated exclusion for SonoVue/Lazarus production domains",
        ),
        "PAIR_REGIME_MAP_PASS": dict(
            status="PASS",
            scopes=["uncoated_sampled_CW"],
            evidence=[evidence(pair)],
            reason="Matched two-radius100/500Pa sampled map only; not50-100kPa coated clusters",
        ),
        "RESPONSE_TABLE_PASS": dict(
            status="PASS" if t.get("status") == "PASS_EMPIRICAL_DOMAIN" else "PARTIAL",
            scopes=["focused_table_numerical"] if t.get("status") == "PASS_EMPIRICAL_DOMAIN" else [],
            evidence=[evidence(table)] if table.exists() else [],
            reason="Numerical interpolation scope is not experimental calibration or a cylinder-wide domain",
        ),
        "TRANSLATION_PASS": dict(
            status="PARTIAL",
            scopes=["solid_shell_sampled_nodes"],
            evidence=[evidence(trans)] if trans.exists() else [],
            reason="Measured 20-cycle solid-shell errors cannot certify SonoVue20ms transport",
        ),
        "LAZARUS_MODEL_SCOPE_PASS": dict(
            status="UNVALIDATED",
            scopes=[],
            evidence=[],
            reason="Independent paper model/rearrangement requires its own benchmark; not a SonoVue table substitution",
        ),
        "LAZARUS_PASS": dict(status="NOT_RUN_GATE", scopes=[], evidence=[]),
        "CYLINDER_PASS": dict(status="NOT_RUN_GATE", scopes=[], evidence=[]),
    }
    registry = root / "audit/PHYSICS_GATES.json"
    atomic_json(
        registry,
        dict(gates=gates, provenance=provenance({}, ["Stage2.5 gates G1-G3"]), zhang_fig3=z["status"]),
    )
    manager = ValidationGateManager(registry)
    for command, folder, scope in [
        ("run-lazarus", "lazarus_paper", "lazarus_paper"),
        ("run-sonovue-clustering", "lazarus_sonovue", "sonovue"),
        ("run-cylinder", "cylinder", "sonovue_cylinder"),
        ("run-frequency-sweep", "frequency_sweep", "sonovue_sweep"),
    ]:
        try:
            manager.check(command, scope=scope)
        except ValidationGateError as error:
            atomic_json(
                root / folder / "STATE.json",
                dict(
                    status="NOT_RUN_GATE",
                    command=command,
                    scope=scope,
                    reason=str(error),
                    UNVALIDATED_OVERRIDE=False,
                ),
            )
    print("Gate decisions saved; no downstream simulation launched", flush=True)


if __name__ == "__main__":
    main()
