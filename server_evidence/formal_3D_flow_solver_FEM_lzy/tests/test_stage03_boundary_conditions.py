from stage03_helpers import *
from fem3d.vascular_diagnostics import PORT_TAGS

def test_named_tags_and_three_natural_outlets_are_used():
    assert PORT_TAGS=={'inlet':4,'outlet_01':3,'outlet_02':5,'outlet_03':2}
    c=config();bc=preflight()['boundary_conditions']
    assert c['boundary_model']['wall']=='u=0'
    assert c['boundary_model']['inlet']=='integral(u dot outward n)=-Q; sigma n=-lambda n'
    assert c['boundary_model']['outlets']=={n:'sigma n=0' for n in PORT_TAGS if n!='inlet'}
    assert bc['essential_boundary'].startswith('WALL only')
