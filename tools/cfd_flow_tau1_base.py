"""Prepare, audit, and resume the research-only fresh tau=1 Base run."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.tau1_base import (  # noqa: E402
    audit_continuous_q_referee,
    audit_full_timestep_replay_8step,
    audit_latest_base_window,
    forensic_dense_discrete_failure,
    forensic_boundary_window_audit,
    prepare_tau1_base,
    run_dense_discrete_diagnostic,
    run_fresh_base,
    salvage_incomplete_dense_diagnostic,
    write_full_timestep_mass_contract,
)
from utils.cfd_flow.full_timestep_evidence import (  # noqa: E402
    write_final_referee_and_acceptance,
    write_instrumented_phase_oracle,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action",
        choices=(
            "prepare",
            "audit-referee",
            "audit-base",
            "forensic",
            "dense-diagnostic",
            "forensic-dense-failure",
            "full-timestep-contract",
            "full-timestep-replay",
            "instrumented-phase-oracle",
            "finalize-full-timestep-referee",
            "salvage-dense",
            "run-base",
        ),
    )
    args = parser.parse_args()
    if args.action == "prepare":
        result = prepare_tau1_base(PROJECT_ROOT)
    elif args.action == "audit-referee":
        result = audit_continuous_q_referee(PROJECT_ROOT)
    elif args.action == "audit-base":
        result = audit_latest_base_window(PROJECT_ROOT)
    elif args.action == "forensic":
        result = forensic_boundary_window_audit(PROJECT_ROOT)
    elif args.action == "dense-diagnostic":
        result = run_dense_discrete_diagnostic(PROJECT_ROOT)
    elif args.action == "forensic-dense-failure":
        result = forensic_dense_discrete_failure(PROJECT_ROOT)
    elif args.action == "full-timestep-contract":
        result = write_full_timestep_mass_contract(PROJECT_ROOT)
    elif args.action == "full-timestep-replay":
        result = audit_full_timestep_replay_8step(PROJECT_ROOT)
    elif args.action == "instrumented-phase-oracle":
        result = write_instrumented_phase_oracle(PROJECT_ROOT)
    elif args.action == "finalize-full-timestep-referee":
        result = write_final_referee_and_acceptance(PROJECT_ROOT)
    elif args.action == "salvage-dense":
        result = salvage_incomplete_dense_diagnostic(PROJECT_ROOT)
    else:
        result = run_fresh_base(PROJECT_ROOT)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("status") not in {"FAIL"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
