"""Frozen pressure design handoff; these are numerical boundary designs, not measured physiology."""
from .network_models.roi_boundary_design import load_frozen_design
from .network_models.parameterized_fem_handoff import measure_extension_geometry, transfer_prediction, export_frozen_design
from .network_models.fem_h0_case import serialize_fixed_pressure_xml, audit_xml_change, create_case
from .paths import project_path

def roi_only_design():
    return load_frozen_design(project_path('results/roi_only/network/roi_boundary_design_frozen.yaml'))
