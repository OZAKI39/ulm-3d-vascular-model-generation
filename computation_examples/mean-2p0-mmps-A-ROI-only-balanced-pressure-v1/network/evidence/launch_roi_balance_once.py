"""Dispatch the unchanged H0 GPU runner exactly once for the pinned ROI-only design case.

Only orchestration/provenance; no numerical options or solver method changes.
Run under a private supervisor with autorestart=false.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def claim_once(path,payload):
    data=json.dumps(payload,indent=2,allow_nan=False)+'\n'
    fd=os.open(path,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o444)
    with os.fdopen(fd,'w') as f:f.write(data);f.flush();os.fsync(f.fileno())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('case','source-case','runner','reference-root'):p.add_argument('--'+arg,type=Path,required=True)
    a=p.parse_args();case=a.case
    pre=json.loads((case/'final_preflight.json').read_text())
    assert pre['status']=='PASS' and pre['planned_CFD_run_count']==1
    assert case.name==pre['case_name']=='mean-2p0-mmps-A-ROI-only-balanced-pressure-v1'
    assert pre['workflow_kind']=='ROI_BOUNDARY_DESIGN' and pre['full_A_network_used_for_boundary_generation'] is False
    assert pre['regression_status']=='PASS' and pre['CFD_feedback_used'] is False
    assert sha(a.runner)==pre['runner_sha256']
    assert sha(case/'roi_boundary_design_frozen.yaml')==pre['frozen_design_sha256']
    assert sha(case/'roi_fem_handoff.json')==pre['handoff_sha256']
    assert json.loads((case/'input_hashes.json').read_text())==pre['input_hashes']
    for name,digest in pre['input_hashes'].items():assert sha(case/name)==digest,name
    diff=pre['configuration_diff']
    assert sha(a.source_case/'run/solver.xml')==diff['original_solver_sha256']
    for name,digest in diff['unchanged_input_hashes'].items():assert sha(a.source_case/name)==sha(case/name)==digest,name
    old,new=ET.parse(a.source_case/'run/solver.xml').getroot(),ET.parse(case/'run/solver.xml').getroot()
    for name,value in pre['applied_cap_pressure_pa'].items():
        path=f".//Add_BC[@name='OUTLET_0{name[1]}']/Value"
        assert float(new.find(path).text)==value
        new.find(path).text=old.find(path).text
    assert ET.tostring(old)==ET.tostring(new),'Non-outlet configuration changed'
    assert not (case/'run/solver.log').exists() and not (case/'run/STOP_SIM').exists()
    command=[sys.executable,'-B',str(a.runner),'--case',str(case),'--reference-root',str(a.reference_root)]
    claim_once(case/'reports/one_scientific_case_dispatch.json',dict(
        status='DISPATCHED',scientific_case_count=1,dispatch_count=1,command=command,
        preflight_sha256=sha(case/'final_preflight.json'),start_unix=time.time(),
        CFD_retuning_permitted=False,automatic_restart_permitted=False))
    code=subprocess.run(command,check=False).returncode
    (case/'reports/dispatch_completion.json').write_text(json.dumps(dict(
        runner_exit_code=code,end_unix=time.time(),new_scientific_CFD_run_count=1,
        scientific_case_count=1,parameter_retuning=False),indent=2)+'\n')
    return code


if __name__=='__main__':raise SystemExit(main())
