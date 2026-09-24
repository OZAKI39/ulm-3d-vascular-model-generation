from copy import deepcopy
import pytest
import yaml
from stage03_helpers import *
from fem3d.vascular import validate_reference

def test_flow_is_read_from_current_source_and_not_experimental():
    c=config();source=yaml.safe_load(Path(c['source']['configuration_path']).read_text())
    close(c['physics']['inlet_volume_flow_m3_s'],source['boundary_conditions']['target_volume_flow_m3_s'])
    assert c['condition_type']=='REFERENCE_NUMERICAL_CONDITION' and c['experimental'] is False
    assert c['warning']=='NOT_EXPERIMENTAL_PUMP_FLOW' and not c['source']['outlet_pressures_imported']
    template=yaml.safe_load((ROOT/'configs/stage03_experiment_template.yaml').read_text())
    assert template['inlet_volume_flow_m3_s'] is None and template['status']=='NOT_EXECUTED'

@pytest.mark.parametrize('change',[{'experimental':True},{'warning':'pump'},{'condition_type':'EXPERIMENTAL'}])
def test_unidentified_or_experimental_condition_is_rejected(change):
    c=deepcopy(config());c.update(change)
    with pytest.raises(ValueError):validate_reference(c)
