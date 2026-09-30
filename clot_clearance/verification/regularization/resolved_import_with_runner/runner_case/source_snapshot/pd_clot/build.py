import json
from pathlib import Path
import numpy as np
from .geometry import make_cloud
from .mechanics import evaluate_cloud, timestep_limit


def main():
    root = Path(__file__).resolve().parents[1]
    c = json.loads((root/'configs/straight_pipe.json').read_text())
    g = make_cloud(c['clot'])
    force, _, J, _, _ = evaluate_cloud(g, g.X.copy(), np.ones(len(g.pairs)), c['material'], c['safety'])
    if np.max(np.abs(force)) > 1e-12 or np.max(np.abs(J-1)) > 1e-10:
        raise RuntimeError('Undeformed state failed')
    print(f'Compiled NOSB-PD: {len(g.X)} particles, {len(g.pairs)} bonds; dt bound {timestep_limit(g,c):.6g} s')


if __name__ == '__main__': main()
