"""Post-audit algebra cross-check using saved states only; no force evaluations."""
from pathlib import Path
import json
import numpy as np
from pd_clot.regularization.equilibrium import make_fixture

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'verification/damage_energy_audit_v1'


def main():
    c = json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text())
    cloud = make_fixture(c)
    i, j = cloud.pairs.T
    rows = []
    for path in sorted((OUT/'actual_run_replay').glob('*/jump_*/particle_audit.npz')):
        z = np.load(path)
        dg = z['D_before'] - z['D_after']
        eta = z['x'][j] - z['x'][i]
        delta = np.zeros_like(z['old_F'])
        for ends, other in [(i,j),(j,i)]:
            residual = eta - np.einsum('bij,bj->bi', z['old_F'][ends], cloud.xi)
            moment = np.einsum('bi,bj,b->bij', residual, cloud.xi, cloud.weight*dg*cloud.volume[other])
            np.add.at(delta, ends, moment)
        predicted = delta @ np.linalg.inv(z['new_K'])
        observed = z['new_F'] - z['old_F']
        rows.append(dict(case=path.parent.parent.name, jump=path.parent.name,
            maximum_absolute_component_error=float(np.max(abs(observed-predicted))),
            maximum_observed_component_change=float(np.max(abs(observed)))))
    result = dict(label='POST-AUDIT ALGEBRA CROSS-CHECK; no new acceptance tolerance',
        identity='F1-F0 = sum_j[w*(g1-g0)*Vj*(eta-F0*xi) outer xi] @ inverse(K1)',
        significance='The weight change acts on the pre-existing non-affine residual; spatial alignment matters. The identity is algebraic and does not itself determine the sign of the energy change.',
        historical_jumps=len(rows),
        maximum_absolute_component_error=max(r['maximum_absolute_component_error'] for r in rows), rows=rows)
    with (OUT/'comparisons/RECONSTRUCTION_IDENTITY.json').open('x') as f:
        json.dump(result,f,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))


if __name__=='__main__':
    main()
