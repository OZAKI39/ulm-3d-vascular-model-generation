from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils.cfd_flow.wall_qvalue_oracle import run_wall_qvalue_oracle  # noqa: E402


if __name__ == "__main__":
    print(json.dumps(run_wall_qvalue_oracle(ROOT), indent=2))
