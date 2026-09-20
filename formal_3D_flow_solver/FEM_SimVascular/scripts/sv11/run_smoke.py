#!/usr/bin/env python3
import difflib, json, shutil, sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv11 import REPORT,load,compare_production_xml
from sv_validation.sv11_runtime import OPTIONS,petsc_ls,run_solver
from sv_validation.provenance import sha256,write_json
from sv_validation.validation import parse_result_vtu

case=ROOT/'outputs/sv1_1/official_fluid_smoke';case.mkdir(parents=True,exist_ok=True)
assert not list(case.glob('4-procs/*.vtu'))
records=[]
for record in load('official_smoke_inputs','sv1')['files']:
    original=ROOT/record['path'];assert sha256(original)==record['sha256']
    name=original.relative_to(ROOT/'outputs/sv1/official_fluid_smoke');target=case/name
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(original,target)
    records.append({'path':str(target.relative_to(ROOT)),'sha256':sha256(target),'source_sha256':record['sha256']})
original=ROOT/'outputs/sv1/official_fluid_smoke/solver.xml'
tree=petsc_ls(ET.parse(original));ET.indent(tree);tree.write(case/'solver.xml',encoding='utf-8',xml_declaration=True)
assert compare_production_xml(original,case/'solver.xml')
write_json(REPORT/'petsc_smoke_inputs.json',{'files':records,'xml_only_LS_changed':True,'options':OPTIONS})
run=run_solver(case,'petsc_smoke',float(tree.findtext('.//Time_step_size')))
files=[];error=None
for path in sorted(case.glob('4-procs/result_*.vtu')):
    try:
        mesh,u,p=parse_result_vtu(path)
        files.append({'path':str(path.relative_to(ROOT)),'sha256':sha256(path),'points':mesh.n_points,'cells':mesh.n_cells,
                      'velocity_finite':bool(np.isfinite(u).all()),'pressure_finite':bool(np.isfinite(p).all()),
                      'velocity_max':float(np.linalg.norm(u,axis=1).max()),'pressure_range':[float(p.min()),float(p.max())]})
    except Exception as exc:error=str(exc)
history=run['history'];reasons=history['petsc_reasons']
passed=run['exit_code']==0 and bool(files) and error is None and not run['monitor_stop'] and bool(reasons) and not any(r['diverged'] for r in reasons) and not history['failed_linear_solves'] and not history['unparsed_rows']
result={'status':'PASS' if passed else 'FAIL','reason':None if passed else 'PETSC_SMOKE_FAIL','run':run,'result_files':files,'field_read_error':error,
        'xml_only_LS_changed':True,'production_run_authorized_by_gate':passed}
write_json(REPORT/'petsc_smoke.json',result);print(json.dumps({k:v for k,v in result.items() if k!='run'},indent=2))
raise SystemExit(0 if passed else 1)
