"""Run the schema-v2 validated Tau1 production CFD entry point."""

from __future__ import annotations

import sys
from pathlib import Path

from utils.cfd_flow.config import load_cfd_flow_config
from utils.cfd_flow.io import FlowError
from utils.cfd_flow.pipeline import SUCCESS_STATUS, print_result, run_cfd_flow


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "cfd_flow.yaml"


def main(argv: list[str] | None = None) -> int:
    """Execute the production flow stage with zero or one YAML argument."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    if len(arguments) > 1:
        print(
            "Usage: python cfd_flow.py "
            "[configs/cfd_flow.yaml|configs/cfd_flow_promotion_regression.yaml]"
        )
        return 1
    config_path = Path(arguments[0]) if arguments else DEFAULT_CONFIG
    try:
        config = load_cfd_flow_config(config_path, project_root=PROJECT_ROOT)
        result = run_cfd_flow(config, project_root=PROJECT_ROOT)
        print_result(result)
        return 0 if result.status == SUCCESS_STATUS else 2
    except (FlowError, FileNotFoundError, ValueError) as error:
        status = error.status if isinstance(error, FlowError) else "CFD_FLOW_INPUT_INVALID"
        print(f"CFD flow failed: {error}")
        print(f"STATUS: {status}")
        print("NEXT: REVIEW CFD FLOW FAILURE EVIDENCE")
        return 1
    except Exception as error:  # keep implementation failures distinguishable
        print(f"CFD flow failed: {error}")
        print("STATUS: CFD_FLOW_INTERNAL_ERROR")
        print("NEXT: REVIEW CFD FLOW INTERNAL ERROR")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
