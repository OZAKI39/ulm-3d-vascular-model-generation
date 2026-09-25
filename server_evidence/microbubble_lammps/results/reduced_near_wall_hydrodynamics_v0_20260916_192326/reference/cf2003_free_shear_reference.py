"""Bounded CF2003 Table17 reference; source-validated fit, not all-gap exactness."""
import json
import math
from pathlib import Path

_data=json.loads(Path(__file__).with_name('CF2003_FREE_SHEAR_REFERENCE_V0.json').read_text())
assert _data['REFERENCE_RESOLUTION_PASS'] and _data['corrigendum_checked']
_domain=tuple(_data['certified_epsilon_domain'])

def _evaluate(epsilon, coefficients):
    e=float(epsilon)
    if not math.isfinite(e) or not _domain[0]<=e<=_domain[1]:
        raise ValueError('CF2003_OUTSIDE_CERTIFIED_DOMAIN')
    x=math.log(e)
    p=coefficients[-1]
    for c in coefficients[-2::-1]: p=p*x+c
    answer=1/p
    if not math.isfinite(answer) or not 0<answer<1:
        raise ArithmeticError('INVALID_CF2003_REFERENCE_VALUE')
    return answer

def CF2003_FU(epsilon):
    return _evaluate(epsilon,_data['u_coefficients'])

def CF2003_FOMEGA(epsilon):
    return _evaluate(epsilon,_data['omega_coefficients'])
