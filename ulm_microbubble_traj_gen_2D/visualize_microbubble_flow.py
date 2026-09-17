"""
Visualize field-based microbubble flow results.
"""

from __future__ import annotations
from pathlib import Path
if __package__ in (None, ""):
    import sys
    sys.path.append(str(Path(__file__).resolve().parents[1]))
    from ulm_microbubble_traj_gen_2D.utils.visualization.results.cli import main, parse_args
else:
    from .utils.visualization.results.cli import main, parse_args


__all__ = ["main", "parse_args"]


if __name__ == "__main__":
    main()
