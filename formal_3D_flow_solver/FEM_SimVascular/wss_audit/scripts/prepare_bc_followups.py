"""Prepare independent one-outlet perturbation inputs; NEVER starts a solver.
The mesh is referenced read-only through a symlink, with content hashes.
Use --source-case /absolute/H0 --output /new/audit/cases on the solver server.
"""
from pathlib import Path
import argparse,json,hashlib,shutil,xml.etree.ElementTree as ET
p=argparse.ArgumentParser();p.add_argument('--source-case',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
base=ET.parse(a.source_case/'run/solver.xml');pp=float(base.find(".//Add_BC[@name='OUTLET_02']/Value").text);p0=float(base.find(".//Add_BC[@name='OUTLET_03']/Value").text)
for sign,name in [(-1,'O2_minus1pct'),(1,'O2_plus1pct')]:
    out=a.output/name
    if out.exists():raise FileExistsError(out)
    (out/'run').mkdir(parents=True);(out/'reports').mkdir();(out/'SV_MESH').symlink_to((a.source_case/'SV_MESH').resolve(),target_is_directory=True)
    for f in ['policy.json','run/PETSC_OPTIONS.txt']:shutil.copy2(a.source_case/f,out/f)
    tr=ET.parse(a.source_case/'run/solver.xml');tr.find(".//Add_BC[@name='OUTLET_02']/Value").text=format(pp+sign*.01*(pp-p0),'.17g');tr.write(out/'run/solver.xml',encoding='utf-8',xml_declaration=True)
    paths=[out/'policy.json',out/'run/PETSC_OPTIONS.txt',out/'run/solver.xml']+sorted((out/'SV_MESH').rglob('*'))
    hashes={str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths if f.is_file()}
    (out/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    (out/'UNRUN.json').write_text(json.dumps(dict(solver_run=False,changed_only='OUTLET_02 pressure',delta_Pa=sign*.01*(pp-p0),reference='1 percent of current O2-O3 pressure difference'),indent=2)+'\n')
    print(out)
