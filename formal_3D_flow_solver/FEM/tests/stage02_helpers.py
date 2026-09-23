import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def read_case(name,kind):
    return json.loads((ROOT/f"outputs/stage02/cases/{name}/qc/{kind}.json").read_text())


def read_report(name):
    return json.loads((ROOT/f"reports/stage02/{name}.json").read_text())


def coefficient_scaling_error(name,factors):
    baseline=np.load(ROOT/"outputs/stage02/cases/pipe_medium_natural/solution/primary_checkpoint.npz")
    actual=np.load(ROOT/f"outputs/stage02/cases/{name}/solution/primary_checkpoint.npz")
    result={}
    for field,factor in zip(("velocity","pressure","lambda_pa"),factors):
        if field=="lambda_pa":
            a,b=actual[field],baseline[field]
        else:
            def ordered(data):
                xyz=data[field+"_coordinates_m"]
                indices=np.lexsort(xyz[:,::-1].T)
                return xyz[indices],data[field+"_values"][indices]
            xyza,a=ordered(actual); xyzb,b=ordered(baseline)
            np.testing.assert_allclose(xyza,xyzb,rtol=0,atol=1e-17)
        result[field]=float(np.linalg.norm(a-factor*b)/np.linalg.norm(factor*b))
    return result
