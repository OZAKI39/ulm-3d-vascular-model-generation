# Ultraliser pipeline

```text
swc_stl_model_generate.py
      |
      v
load and validate configs/swc_stl_model_generate.yaml
      |
      v
load saved ROI / ROIRecord
      |
      v
ultraliser_backend.py
      |
      +--> export canonical SWC (source radii unchanged)
      |
      +--> create vascular H5
      |      feed radius = source radius * 0.91
      |
      +--> call official ultraVessMorpho2Mesh once
      |
      v
ultraliser_qc.py
      |
      +--> export STL um / VTP um / STL m
      |
      +--> topology QC
      |
      +--> source-radius fidelity QC
      |
      v
final output and reconstruction report
```

## Core responsibilities

- `swc_stl_model_generate.py`: YAML-driven saved-ROI selection and orchestration; it does not
  maintain per-parameter command-line overrides.
- `utils/cfd_lumen/model_yaml_config.py`: strict YAML schema, path resolution, selector checks,
  and construction of the validated reconstruction/QC configuration.
- `utils/cfd_lumen/ultraliser_backend.py`: ROI validation and serialization, canonical SWC export,
  H5 adapter and radius mapping, executable discovery, official subprocess invocation, and explicit
  failure reporting.
- `utils/cfd_lumen/ultraliser_qc.py`: source branch extraction, STL/VTP unit conversion, whole-mesh
  topology checks, radius cross-section fidelity, and the concise reconstruction report.
- `utils/cfd_lumen/roi_io.py` and `utils/sampling/`: provide saved `ROIRecord` objects and preserve
  source centerlines, radii, `CUT_PORT` metadata, and local/global mappings.
- `Ultraliser/`: upstream C++ `ultraVessMorpho2Mesh`; its reconstruction algorithm is not modified.

There is no alternate surface backend and no geometry fallback. A failed Ultraliser invocation is
reported as `ULTRALISER_EXECUTION_FAILED`.
