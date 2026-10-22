"""Clone a solved control case and change only O2 by +/-1% of (O2-O3)."""
from pathlib import Path
import argparse,json,shutil,xml.etree.ElementTree as ET
from case_common import *
assert not (V/'stage4/user_cancellation.json').exists(), 'SKIPPED_BY_USER: pressure cases must not be generated'
p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,default=V/'stage3/vessel_medium');args=p.parse_args();base=args.baseline.resolve();assert json.loads((base/'reports/execution.json').read_text())['status']=='PASS';assert json.loads((base/'reports/flow_quality.json').read_text())['accepted_final_and_log_checks']
root=ET.parse(base/'run/solver.xml');p0={n:float(root.find('.//Add_BC[@name="%s"]/Value'%n).text) for n in ['OUTLET_01','OUTLET_02','OUTLET_03']};delta=.01*(p0['OUTLET_02']-p0['OUTLET_03']);assert delta>0
for sign,name in [(-1,'O2_minus1pct'),(1,'O2_plus1pct')]:
 case=V/'stage4'/name;assert not case.exists();case.mkdir(parents=True);shutil.copytree(base/'SV_MESH',case/'SV_MESH');(case/'run').mkdir();shutil.copy2(base/'run/PETSC_OPTIONS.txt',case/'run/PETSC_OPTIONS.txt')
 policy=json.loads((base/'policy.json').read_text());policy.update(case=name,kind='fixed_grid_outlet_pressure_CFD',mesh_reference=base.name,comparison_baseline_case=base.name,perturbation=dict(parameter='OUTLET_02_pressure_Pa',baseline=p0['OUTLET_02'],delta_Pa=sign*delta,fraction_of_O2_minus_O3=sign*.01));dump(case/'policy.json',policy)
 tree=ET.parse(base/'run/solver.xml');tree.find('.//Add_BC[@name="OUTLET_02"]/Value').text=repr(p0['OUTLET_02']+sign*delta);tree.write(case/'run/solver.xml',encoding='utf-8',xml_declaration=True)
 changes=[]
 for a,b in zip(root.getroot().iter(),tree.getroot().iter()):
  assert a.tag==b.tag and a.attrib==b.attrib
  if (a.text or '').strip()!=(b.text or '').strip():changes.append(dict(tag=a.tag,before=a.text,after=b.text))
 assert len(changes)==1 and changes[0]['tag']=='Value';lock_case(case);dump(case/'reports/single_factor_check.json',dict(baseline_case=str(base),semantic_XML_changes=changes,identical_PETSC_OPTIONS=sha(base/'run/PETSC_OPTIONS.txt')==sha(case/'run/PETSC_OPTIONS.txt'),all_mesh_files_identical=all(sha(p)==sha(case/'SV_MESH'/p.relative_to(base/'SV_MESH')) for p in (base/'SV_MESH').rglob('*') if p.is_file())));print(name,changes)
dump(V/'stage4/boundary_group.json',dict(baseline_case=base.name,baseline_path=str(base),mesh=base.name,O1_Pa=p0['OUTLET_01'],O2_Pa=p0['OUTLET_02'],O3_Pa=p0['OUTLET_03'],delta_Pa=delta,dt_s=policy['dt_s'],MPI_ranks=policy['MPI_ranks'],only_factor='O2 minus O3 pressure difference'))
