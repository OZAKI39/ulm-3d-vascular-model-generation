"""Actual entry must reject a self-consistent but altered frozen-case material."""
import sys,json,tempfile,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
V=Path(__file__).resolve().parents[1];C=V.parent;H=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1')
(V/'runtime').mkdir(exist_ok=True)
with tempfile.TemporaryDirectory(dir=V/'runtime',prefix='negative_frozen_config_') as tmp:
 p=Path(tmp);(p/'run').mkdir();tree=ET.parse(H/'run/solver.xml');mu=tree.find('.//Viscosity/Value');mu.text=str(2*float(mu.text));tree.write(p/'run/solver.xml',encoding='utf-8',xml_declaration=True)
 policy=json.loads((H/'policy.json').read_text());policy['nu_m2_s']*=2;(p/'policy.json').write_text(json.dumps(policy));(p/'frozen_flow').symlink_to(H/'frozen_flow',target_is_directory=True);(p/'SV_MESH').symlink_to(H/'SV_MESH',target_is_directory=True)
 r=subprocess.run([sys.executable,'-B',str(C/'rotate_visualization/prepare_surface_data.py'),'--case',str(p),'--output',str(p/'must_not_exist')],capture_output=True,text=True)
 assert r.returncode!=0 and 'Configuration differs from frozen solve' in r.stderr and not (p/'must_not_exist').exists()
 record=dict(test='consistent mu and nu edits on a case with original frozen solution',altered_mu_Pa_s=float(mu.text),altered_nu_m2_s=policy['nu_m2_s'],exit_code=r.returncode,expected_rejection=True,stderr=r.stderr,formal_case_untouched=True)
(V/'stage1/data/frozen_configuration_negative_test.json').write_text(json.dumps(record,indent=2)+'\n');print('Frozen configuration mismatch rejected by actual production entry.')
