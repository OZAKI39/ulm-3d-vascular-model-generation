#!/usr/bin/env python3
"""Freeze the successful one-rank mathematical checks before derived fields."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"tests"))
from fem3d.audit import timestamp,write_json
from stage02_helpers import coefficient_scaling_error

report=ROOT/"reports/stage02"
for filename in ("canonical_core_pytest.xml","scaling_pytest.xml"):
    t=ET.parse(report/filename).getroot().find("testsuite").attrib
    assert int(t["failures"])==int(t["errors"])==int(t["skipped"])==0
scaling={}
for name,factors in {"q_half":(.5,.5,.5),"q_double":(2,2,2),"reverse":(-1,-1,-1),"mu_double":(1,2,2),"rho_double":(1,1,1)}.items():
    scaling[name]=coefficient_scaling_error(name,factors)
    assert max(scaling[name].values())<1e-11
assert max(scaling["rho_double"].values())==0
case_names=[f"{profile}_{mode}" for profile in ("pipe_coarse","pipe_medium","pipe_fine") for mode in ("natural","reference")]+list(scaling)
for name in case_names:
    b=ROOT/"outputs/stage02/cases"/name
    f=json.loads((b/"qc/flux.json").read_text())
    s=json.loads((b/"qc/solver.json").read_text())
    assert f["status"]=="PASS" and s["converged_reason"]>0 and s["physical_fields_finite"]
    assert s["real_global_dofs"]==1 and s["boundary_conditions"]["pressure_dirichlet_dofs"]==0
result={"status":"PASS","timestamp":timestamp(),"case_names":case_names,"case_count":len(case_names),
    "scaling_relative_errors":scaling,"formal_boundary_preserved":True,
    "analytic_reference":"Additional traction-compatible case explicitly authorized by user",
    "authorization_for_derived_export":"All one-rank core mathematical tests passed before postprocessing implementation",
    "stage3_started":False}
write_json(report/"core_validation.json",result)
write_json(ROOT/"inputs/stage02/core_gate.json",result)
print(json.dumps(result,indent=2))
