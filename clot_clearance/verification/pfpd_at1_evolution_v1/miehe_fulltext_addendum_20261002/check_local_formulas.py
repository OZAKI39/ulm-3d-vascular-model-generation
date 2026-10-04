"""Document-only scalar consistency checks of the supplied Miehe PDF.

No project solver, saved simulation, load schedule or material parameter is
modified. The unit coefficient A=1 below is an algebraic normalization,
not a new Gc or a calibrated material parameter.
"""
from pathlib import Path
import json
from scipy.optimize import minimize_scalar


def minimum(objective):
    result = minimize_scalar(objective, bounds=(0., 1.), method='bounded',
                             options={'xatol': 1e-14})
    candidates = [(0., objective(0.)), (float(result.x), float(result.fun)), (1., objective(1.))]
    return dict(raw_s=float(result.x), success=bool(result.success),
                selected_s=min(candidates, key=lambda pair: pair[1])[0], candidates=candidates)


def main():
    rows = []
    for ratio in [0., .01, .5, 1., 2., 10.]:
        # A = gc/(2*ell); ratio = H/A; eta=0 and spatially uniform phase.
        d_balance = ratio/(1+ratio)
        d_at1 = 0. if ratio <= 1 else 1-1/ratio
        plus = minimum(lambda d: d*d+(1-d)**2*ratio)
        printed_minus = minimum(lambda d: d*d-(1-d)**2*ratio)
        at1 = minimum(lambda d: 2*d+(1-d)**2*ratio)
        rows.append(dict(H_over_A=ratio, Miehe_Eq41_local_solution=d_balance,
            positive_history_energy_minimum=plus,
            literal_printed_Eq55_energy_minimum=printed_minus,
            positive_energy_derivative_at_Eq41_solution=2*d_balance-2*(1-d_balance)*ratio,
            printed_minus_derivative_at_Eq41_solution=2*d_balance+2*(1-d_balance)*ratio,
            AT1_analytic_for_same_numerical_coefficient=d_at1,
            AT1_independent_energy_minimum=at1,
            algebraic_comparison_only_not_identification_of_physical_normalizations=True))
    output = dict(scope='Scalar algebra audit only; no mechanics or coupled rerun',
        source=dict(pdf='1-s2.0-S0045782510001283-main.pdf',
                    doi='10.1016/j.cma.2010.04.011',
                    sha256='9fa6625138ee7a5ed073f5f7bbeaae887f078a4a4be848e3e9a8771f3b77bdb1'),
        dimensionless_normalization='A=gc/(2*ell); divide local densities by A',
        comparison_formulas=dict(Miehe_Eq41='d = H/(H + gc/(2*ell))',
            balance_consistent_energy='(1-d)^2 H + gc/(2*ell) d^2',
            literal_printed_Eq55='-(1-d)^2 H + gc/(2*ell) d^2 [eta=0, gradient=0]',
            current_AT1='(1-s)^2 Y + 2 Yc s',
            hypothetical_quadratic_at1_same_prefactor='(1-s)^2 Y + 2 Yc s^2 -> s=Y/(Y+2Yc)'),
        publisher_correction_confirmed=False, rows=rows)
    Path(__file__).with_name('LOCAL_FORMULA_CHECK.json').write_text(json.dumps(output,indent=2)+'\n')
    assert all(r['positive_history_energy_minimum']['success'] and
               r['literal_printed_Eq55_energy_minimum']['success'] and
               r['AT1_independent_energy_minimum']['success'] for r in rows)
    assert max(abs(r['positive_history_energy_minimum']['selected_s']-
                   r['Miehe_Eq41_local_solution']) for r in rows)<5e-8
    assert all(r['literal_printed_Eq55_energy_minimum']['selected_s']==0 for r in rows)
    assert max(abs(r['AT1_independent_energy_minimum']['selected_s']-
                   r['AT1_analytic_for_same_numerical_coefficient']) for r in rows)<5e-8
    print('Six normalized scalar cases checked. Printed minus sign does not reproduce Eq.41.')


if __name__ == '__main__': main()
