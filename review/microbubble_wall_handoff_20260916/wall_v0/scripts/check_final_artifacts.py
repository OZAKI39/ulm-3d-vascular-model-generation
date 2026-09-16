from pathlib import Path
import json,hashlib,re,subprocess
from PIL import Image
import vtk
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for p in (R/'src').glob('*.py'):compile(p.read_text(),str(p),'exec')
for p in (R/'scripts').glob('*.py'):compile(p.read_text(),str(p),'exec')
v=json.loads((R/'VISUALIZATION_PROVENANCE.json').read_text());assert v['human_visual_review']=='PENDING';assert sha(R/'scripts/generate_wall_visualizations.py')==v['script_sha256']
for path,h in {**v['sources'],**v['outputs']}.items():assert sha(R/path)==h,path
images=list((R/'visualization').glob('VIS_*.png'));assert len(images)==12
for p in images:
 with Image.open(p) as im:assert im.width>=800 and im.height>=400;im.verify()
counts={};required=['particle_id','time_s','diameter_um','velocity_m_s','angular_velocity_rad_s','gap_m','epsilon','wall_normal','triangle_id','rms_over_a','normal_spread_deg','wall_active','constraint_active']
for p in (R/'visualization/paraview').glob('*.vtp'):
 reader=vtk.vtkXMLPolyDataReader();reader.SetFileName(str(p));reader.Update();poly=reader.GetOutput();counts[p.name]=poly.GetNumberOfPoints()
 if p.name!='REAL_SURFACE_VALIDITY_SAMPLES.vtp':
  names=[poly.GetPointData().GetArrayName(i) for i in range(poly.GetPointData().GetNumberOfArrays())];assert all(k in names for k in required),(p,names)
assert counts['CASE_I_TRAJECTORY.vtp']==0 and counts['CASE_J_TRAJECTORIES.vtp']==8 and counts['LOCAL_PLANE_VALIDITY.vtp']==10000
for p in [R/'README.md',R/'MICROBUBBLE_WALL_HYDRODYNAMICS_V0_REPORT.md']:
 for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if not target.startswith('http'):assert (p.parent/target).exists(),target
assert json.loads((R/'provenance/REMOTE_CONTROL_CLOSED.json').read_text())['status']=='PASS'
assert json.loads((R/'FINAL_STATUS.json').read_text())['MICROBUBBLE_WALL_HYDRODYNAMICS_V0']=='BLOCKED_LOCAL_PLANE_VALIDITY'
p=subprocess.run(['sha256sum','-c','--quiet','REMOTE_SHA256SUMS'],cwd=R,capture_output=True,text=True);assert p.returncode==0,(p.stdout,p.stderr)
result={'status':'PASS','PNG_count':12,'VTK_points':counts,'visualization_sources_and_outputs_SHA':'PASS','python_syntax':'PASS','relative_document_links':'PASS','remote_download_manifest_recheck':'PASS','human_visual_review':'PENDING'};(R/'validation/FINAL_ARTIFACT_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
