"""Continue adaptive-flux steady CFD once, then run the exact audit on PASS."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.adaptive_flux_steady import STEADY_PENDING_AUDIT  # noqa: E402
from utils.cfd_flow.adaptive_flux_steady_exact_audit import (  # noqa: E402
    run_adaptive_flux_steady_exact_audit,
)
from utils.cfd_flow.adaptive_flux_steady_resume import (  # noqa: E402
    STEADY_PASS_0P5,
    finalize_resume_with_exact_audit,
    run_adaptive_flux_steady_resume,
)
from utils.cfd_flow.io import write_json  # noqa: E402


def main() -> int:
    steady = run_adaptive_flux_steady_resume(PROJECT_ROOT)
    exact = None
    if steady.get("status") == STEADY_PENDING_AUDIT:
        exact = run_adaptive_flux_steady_exact_audit(
            PROJECT_ROOT, Path(steady["run_root"])
        )
    result = finalize_resume_with_exact_audit(steady=steady, exact_audit=exact)
    write_json(
        Path(steady["run_root"]) / "qc" / "adaptive_flux_resume_0p5_final_manifest.json",
        result,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("status") == STEADY_PASS_0P5 else 1


if __name__ == "__main__":
    raise SystemExit(main())
