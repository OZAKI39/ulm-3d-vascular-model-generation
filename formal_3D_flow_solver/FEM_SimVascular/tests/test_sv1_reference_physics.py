import yaml
from sv_validation.validation import reference_condition

def test_frozen_si_numerical_reference(root, report):
    cfg=yaml.safe_load((root/"configs/sv_reference.yaml").read_text())
    p=reference_condition(cfg)
    source=yaml.safe_load((root/"inputs/fem_reference/stage03_reference_vascular.yaml").read_text())
    assert cfg == report("reference_physics")
    assert p["density_kg_m3"] == source["physics"]["density_kg_m3"]
    assert p["inlet_volume_flow_m3_s"] == source["physics"]["inlet_volume_flow_m3_s"]
