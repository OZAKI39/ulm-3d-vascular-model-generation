import xml.etree.ElementTree as ET
import pytest
from sv_validation.sv11 import ROOT,REPORT,compare_production_xml,validate_petsc_xml
from sv_validation.validation import ValidationError

def candidate():
    path=ROOT/'configs/sv1_1/sv_flow_petsc.xml'
    if not path.exists():pytest.skip('Production XML forbidden before smoke passes')
    return path

def test_production_only_linear_algebra_changes():
    assert compare_production_xml(ROOT/'configs/sv_flow.xml',candidate())
    assert validate_petsc_xml(candidate())

@pytest.mark.parametrize('xpath',['.//Mesh_file_path','.//Face_file_path','.//Time_step_size','.//Density','.//Viscosity/Value','.//Add_BC/Value','.//Add_BC/Profile','.//Increment_in_saving_VTK_files'])
def test_frozen_production_changes_rejected(tmp_path,xpath):
    tree=ET.parse(candidate());tree.find(xpath).text='changed'
    path=tmp_path/'changed.xml';tree.write(path)
    with pytest.raises(ValidationError):compare_production_xml(ROOT/'configs/sv_flow.xml',path)

@pytest.mark.parametrize('tag',['Tolerance','Absolute_tolerance','Krylov_space_dimension','PETSc_options','Fieldsplit'])
def test_ignored_xml_options_rejected(tmp_path,tag):
    tree=ET.parse(candidate());ET.SubElement(tree.find('.//LS'),tag).text='claimed enabled'
    path=tmp_path/'ignored.xml';tree.write(path)
    with pytest.raises(ValidationError,match='Ignored'):validate_petsc_xml(path)
