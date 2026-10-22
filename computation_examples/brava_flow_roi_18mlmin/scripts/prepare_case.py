from pathlib import Path
import argparse,json,hashlib,xml.etree.ElementTree as ET,shutil
ROOT=Path(__file__).resolve().parents[1]
TEMPLATE=Path('/home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1')
p=argparse.ArgumentParser();p.add_argument('name');p.add_argument('--mode',choices=['calibration','pressure'],required=True);p.add_argument('--pressures',type=Path);args=p.parse_args()
case=ROOT/'cases'/args.name;assert not case.exists();(case/'run').mkdir(parents=True);(case/'reports').mkdir()
(case/'SV_MESH').symlink_to('../../mesh/SV_MESH',target_is_directory=True)
geom=json.loads((ROOT/'reports/geometry_qc.json').read_text());g=json.loads((ROOT/'inputs/geometry.json').read_text())
Q=3e-7;area=geom['ports']['INLET']['actual']['area_m2'];U=Q/area;dt=.01
policy=dict(dt_s=dt,save_interval_steps=10,steady_last_intervals=5,velocity_change_limit=1e-5,flow_change_limit=1e-6,mass_limit=1e-6,
 maximum_total_steps=500,maximum_wall_time_s=21600,Q_target_m3_s=Q,Umean_m_s=U,A_in_m2=area,
 initial_state='zero',dt_rationale='Implicit generalized-alpha march toward a stationary solution; 0.01 s is a steady relaxation time step, not a time-resolved physiological waveform claim. Strict consecutive steady and residual gates retained.',
 rho_kg_m3=1056,mu_Pa_s=.00345312,mode=args.mode)
pressures=json.loads(args.pressures.read_text()) if args.pressures else {'OUTLET_01':0.,'OUTLET_02':0.,'OUTLET_03':0.}
xml=ET.parse(TEMPLATE/'run/solver.xml');xml.find('.//Time_step_size').text=str(dt);xml.find('.//Number_of_time_steps').text='500'
for bc in xml.findall('.//Add_BC'):
 role=bc.get('name');value=bc.find('Value')
 if role=='INLET':value.text=format(-Q,'.17g')
 elif role.startswith('OUTLET_'):
  if args.mode=='calibration' and role in ('OUTLET_01','OUTLET_03'):
   bc.find('Type').text='Dir';value.text='1e-7'
   for tag,text in [('Profile','Parabolic'),('Impose_flux','true'),('Zero_out_perimeter','true')]:ET.SubElement(bc,tag).text=text
  else:value.text=format(pressures[role],'.17g')
ET.indent(xml,space='  ');xml.write(case/'run/solver.xml',encoding='utf-8',xml_declaration=True)
(case/'policy.json').write_text(json.dumps(policy,indent=2)+'\n')
shutil.copy2(TEMPLATE/'run/PETSC_OPTIONS.txt',case/'run/PETSC_OPTIONS.txt')
hashes={str(f.relative_to(case)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [case/'run/solver.xml',case/'run/PETSC_OPTIONS.txt',case/'policy.json']}
for f in (ROOT/'mesh/SV_MESH').rglob('*'):
 if f.is_file():hashes['SV_MESH/'+str(f.relative_to(ROOT/'mesh/SV_MESH'))]=hashlib.sha256(f.read_bytes()).hexdigest()
(case/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
(case/'purpose.json').write_text(json.dumps(dict(mode=args.mode,
 calibration='Qin=3e-7, QO1=QO3=1e-7 m3/s, PO2=0 Pa. Use measured terminal pressures as initial inverse design for a separate final case with natural pressure conditions at all three outlets.' if args.mode=='calibration' else 'Final all-pressure-outlet forward verification; outlet flow rates are measured, not prescribed.',
 prescribed_pressures_Pa=pressures if args.mode=='pressure' else {'OUTLET_02':0},source='real BraVa accepted four-port core'),indent=2)+'\n')
print(case)
