import xml.etree.ElementTree as ET

def test_actual_three_zero_natural_outlets(root):
    equation=ET.parse(root/"configs/sv_flow.xml").find(".//Add_equation")
    bc={b.attrib["name"]: b for b in equation.findall("Add_BC")}
    assert set(bc) == {"WALL","INLET","OUTLET_01","OUTLET_02","OUTLET_03"}
    for role in ("OUTLET_01","OUTLET_02","OUTLET_03"):
        assert bc[role].findtext("Type") == "Neu"
        assert float(bc[role].findtext("Value")) == 0
    assert bc["INLET"].findtext("Profile") == "Flat"
    assert bc["WALL"].findtext("Type") == "Dir" and float(bc["WALL"].findtext("Value")) == 0
