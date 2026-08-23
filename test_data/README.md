# Minimal test-data fixture

This directory contains only the data needed for the active ROI 003274 reconstruction test.

| File group | Rows/items | Purpose |
|---|---:|---|
| `sampling/roi_library/*.npz` | 1 ROI | Local centerline, radius, topology, and CUT_PORT arrays |
| `candidate_rois.csv` | 1 ROI | Loader-compatible candidate manifest |
| `selected_rois.csv` | 1 ROI | Loader-compatible selected manifest |
| `cut_ports.csv` | 3 ports | Auditable ROI boundary intersections |
| `roi_features.csv` | 1 ROI | Saved radius and structural descriptors |
| `global_edges.csv` | 7,418 edges | Exact global edge identity verification for the one source model |
| `analysis_swc_single_component.npz` | 1 component | Global context extension source (7,419 nodes) |

Not included: raw TIFF images, segmentation masks, the original full SWC collection, other ROI
candidates, any other mouse-brain tiles, or generated lumen meshes/reports.

`SHA256SUMS` records the reviewed fixture hashes.
