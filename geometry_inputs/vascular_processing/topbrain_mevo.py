"""Exact MCA M2/M3 masks; other labels remain in the generic registry."""
from dataclasses import dataclass, field
import numpy as np

from .topbrain_qc import TopBrainError, components, contact, physical_points

EXACT = "EXACT_TOPBRAIN_LABEL"
DERIVED = "DERIVED_FROM_TOPBRAIN_LABELS"
COARSE = "COARSE_TOPBRAIN_LABEL"


@dataclass
class TopBrainROI:
    name: str
    mask: np.ndarray | None
    side: str
    territory: str
    segments: list[str]
    native_labels: list[str]
    anatomical_status: str
    status: str = "PASS"
    warnings: list[str] = field(default_factory=list)
    qc: dict = field(default_factory=dict)
    graph: object = None
    modelable_graph: object = None
    output_dir: object = None


def extract_rois(case):
    rois = []
    for side in ("R","L"):
        names = [f"{side}-M2",f"{side}-M3"]
        roi = TopBrainROI(f"{side}MCA_M2M3",None,side,"MCA",["M2","M3"],names,EXACT)
        try:
            masks = [case.mask(n) for n in names]
            roi.mask = masks[0] | masks[1]
            _,roi.qc = components(roi.mask,case.spacing_mm)
            roi.qc["native_voxel_counts"] = {n:int(m.sum()) for n,m in zip(names,masks)}
            if not all(m.any() for m in masks):
                raise TopBrainError("EMPTY_ROI" if not roi.mask.any() else "MISSING_SEGMENT",
                                    "Both native M2 and M3 must be present for complete M2+M3; available union is preserved")
            roi.qc["m2_m3_contact_voxels"] = int(contact(*masks).sum())
            if roi.qc["component_count"] != 1:
                roi.status = "MULTIPLE_COMPONENTS_NEEDS_REVIEW"
                roi.warnings.append("MULTIPLE_COMPONENTS: strict mask retains every component")
            if not roi.qc["m2_m3_contact_voxels"]:
                roi.status = "TOPOLOGY_FRAGMENTED"
                roi.warnings.append("TOPOLOGY_FRAGMENTED: M2 and M3 do not touch")
        except TopBrainError as exc:
            roi.status = exc.code
            roi.warnings.append(str(exc))
        rois.append(roi)
    return rois
