from copy import deepcopy
import pytest
import yaml
from stage03_helpers import *
from fem3d.vascular import validate_reference

def test_mu_uses_source_density_times_kinematic_viscosity():
    c=config();p=c['physics'];source=yaml.safe_load(Path(c['source']['configuration_path']).read_text())['physics']
    close(p['density_kg_m3'],source['density_kg_m3']);close(p['kinematic_viscosity_m2_s'],source['kinematic_viscosity_m2_s'])
    close(p['dynamic_viscosity_pa_s'],p['density_kg_m3']*p['kinematic_viscosity_m2_s'])
    c=deepcopy(c);c['physics']['dynamic_viscosity_pa_s']=c['physics']['kinematic_viscosity_m2_s']
    with pytest.raises(ValueError):validate_reference(c)
