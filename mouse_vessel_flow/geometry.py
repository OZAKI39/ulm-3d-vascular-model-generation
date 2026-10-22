"""Mouse SWC -> connected ROI -> lumen; SWC coordinates/radii are micrometres before SI export.

The retained flow ROI is anchor 003274 of fMOST_0_5_6_0_0_6_0001_02_01.
Ultraliser uses 0.91*source radius (H5 fourth column is diameter); the 1D resistance
network retains original source radii. These conventions must not be conflated.
"""
from .geometry_tools.rodent_vasculature import run_rodent_vasculature_pipeline
from .geometry_tools.sampling.pipeline import run_sampling_from_rodent_run
from .geometry_tools.cfd_lumen import load_sampling_rois
from .geometry_tools.cfd_lumen.ultraliser_backend import run_ultraliser_reconstruction
from .network_models.swc_graph import load_swc, audit_graph
from .paths import project_path

def retained_roi():
    rois=load_sampling_rois(project_path('inputs/geometry/outputs/sampling/20260825_133201_radius_plus_structure_k5'),selected_only=True)
    selected=[r for r in rois if int(r.anchor_id)==3274]
    if len(selected)!=1:raise ValueError('Retained flow ROI identity is ambiguous')
    return selected[0]
