import xml.etree.ElementTree as ET
import pytest
from sv_validation.sv12 import ROOT,CONFIG,OUTPUT,validate_frozen,validate_xml,load
from sv12_support import artifact
import sv_validation.sv12 as sv12
from sv_validation.validation import ValidationError

def test_actual_all_frozen_inputs():assert validate_frozen() and validate_xml(CONFIG/'sv_flow.xml')

def test_actual_running_case_xml_has_only_authorized_differences():
    flow=artifact('flow_execution')
    assert validate_xml(OUTPUT/'vascular_flow/solver.xml',extension=bool(flow['extension_used']))

def test_changed_mesh_hash_rejected(tmp_path,monkeypatch):
    record=next(r for r in load('frozen_input_manifest')['files'] if r['path'].endswith('mesh-complete.mesh.vtu'))
    p=tmp_path/'changed.vtu';p.write_bytes(b'changed mesh')
    manifest={'files':[{**record,'path':str(p)}]}
    monkeypatch.setattr(sv12,'load',lambda name,*args:manifest if name=='frozen_input_manifest' else load(name,*args))
    with pytest.raises(ValidationError,match='FROZEN_INPUT_CHANGED'):validate_frozen()

@pytest.mark.parametrize('tag',['Time_step_size','Density','Add_BC/Value','Viscosity/Value','Mesh_file_path'])
def test_frozen_model_change_rejected(tmp_path,tag):
    tree=ET.parse(CONFIG/'sv_flow.xml');tree.find('.//'+tag).text='changed';p=tmp_path/'changed.xml';tree.write(p)
    with pytest.raises(ValidationError,match='FROZEN_INPUT_CHANGED'):validate_xml(p)
