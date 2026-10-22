"""Fail-closed BraVa calibration -> measured pressure boundary -> forward solve."""
from pathlib import Path
import sys,json,time,subprocess,shutil,hashlib,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];PY=sys.executable

def run(*args):subprocess.run([PY,'-B',*map(str,args)],check=True,cwd=ROOT)
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
while not (ROOT/'cases/diagnostic_gpu8_host_sync/reports/execution.json').exists():time.sleep(5)
run(ROOT/'scripts/validate_gpu_backend.py')
cal=ROOT/'cases/calibration_gpu8_production'
run(ROOT/'scripts/solve_remote.py','--case',cal,'--reference-root',ROOT)
run(ROOT/'scripts/postprocess_flow.py',cal)
pressures=json.loads((cal/'postprocess/designed_outlet_pressures_Pa.json').read_text())
final=ROOT/'cases/balanced_pressure_final';assert not final.exists()
(final/'run').mkdir(parents=True);(final/'reports').mkdir();(final/'SV_MESH').symlink_to('../../mesh/SV_MESH',target_is_directory=True)
xml=ET.parse(cal/'run/solver.xml')
for bc in xml.findall('.//Add_BC'):
 role=bc.get('name')
 if role.startswith('OUTLET_'):
  bc.find('Type').text='Neu';bc.find('Value').text=format(pressures[role],'.17g')
  for tag in ['Profile','Impose_flux','Zero_out_perimeter']:
   element=bc.find(tag)
   if element is not None:bc.remove(element)
ET.indent(xml,space='  ');xml.write(final/'run/solver.xml',encoding='utf-8',xml_declaration=True)
shutil.copy2(cal/'run/PETSC_OPTIONS.txt',final/'run/PETSC_OPTIONS.txt')
policy=json.loads((cal/'policy.json').read_text());policy['mode']='pressure';dump(final/'policy.json',policy)
dump(final/'purpose.json',{'mode':'pressure','boundary_pressures_Pa':pressures,'source':str(cal/'postprocess/designed_outlet_pressures_Pa.json'),'target_split_only':'Flow at every outlet is measured in the forward CFD, never prescribed','initial_state':'zero'})
files=[final/'run/solver.xml',final/'run/PETSC_OPTIONS.txt',final/'policy.json']
hashes={str(p.relative_to(final)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
for p in (ROOT/'mesh/SV_MESH').rglob('*'):
 if p.is_file():hashes['SV_MESH/'+str(p.relative_to(ROOT/'mesh/SV_MESH'))]=hashlib.sha256(p.read_bytes()).hexdigest()
dump(final/'input_hashes.json',hashes)
run(ROOT/'scripts/solve_remote.py','--case',final,'--reference-root',ROOT)
run(ROOT/'scripts/postprocess_flow.py',final)
dump(ROOT/'reports/FLOW_SEQUENCE_COMPLETE.json',{'status':'COMPLETE','case':str(final),'subsequent_steps':'Flow animation, then microbubble calculation; no automatic use of failed output'})
