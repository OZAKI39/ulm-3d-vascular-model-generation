from copy import deepcopy
from stage03_helpers import *
from fem3d.vascular import reynolds

def test_re_uses_projected_area_perimeter_and_actual_reference_properties():
    c=config();r=reynolds(c);g=c['inlet_geometry'];p=c['physics']
    close(r['hydraulic_diameter_m'],4*g['projected_area_m2']/g['projected_perimeter_m'])
    close(r['Re'],p['density_kg_m3']*p['inlet_volume_flow_m3_s']/g['projected_area_m2']*r['hydraulic_diameter_m']/p['dynamic_viscosity_pa_s'])
    changed=deepcopy(c);changed['physics']['inlet_volume_flow_m3_s']*=2/r['Re']
    assert reynolds(changed)['status']=='STOKES_PHYSICS_REVIEW_REQUIRED'
