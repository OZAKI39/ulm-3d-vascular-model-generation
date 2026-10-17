# Methods, conventions and primary sources

Image x points right and y down. Coordinates are ROI-local pixels. The fixed geometric ROI center xc is the reference for v0, and is output separately from the visible candidate mean. No microscope scale is available: v0 in px/s, gradients in 1/s. Positive omega means clockwise rotation on a displayed image.

`v(x) = v0 + A(x-xc)`; `delta=a11+a22`; `omega=a21-a12`; `shear1=(a11-a22)/2`; `shear2=(a12+a21)/2`. Local rigid angular velocity is omega/2. Delta is 2D divergence, not a 3D contact/aggregation rate.

Fit endpoint displacements using SciPy least_squares with Huber loss, normalized spatial design, then form H=[[F,b],[0,0,1]]. Compute L=logm(H)/actual_dt, A=L[:2,:2], v0=L[:2,2]+A*xc. Reject nonpositive time, noninvertible/ill-conditioned or complex-log transforms. This is an interval-equivalent constant generator, not instantaneous motion during acceleration. An analytic mixed finite rotation/translation test verifies the conversion; simply using (F-I)/dt would have biased it.

Require >=8 vectors, rank 3, normalized design condition <=15, hull/ROI area >=0.18, RMS residual <=1.5 px. Reject low texture, peak quality and local median outliers before fitting. Four spatial quadrants are reported; bootstrap uses quadrant resampling (24 nominal replicates) because overlapping PIV vectors are correlated. Measured coverage is below 95%; these intervals are diagnostic and **not calibrated confidence bounds**. Numerical support checks still admit some physically wrong flow fields.

Candidate means/covariances describe visible support. In general dSigma/dt=A Sigma+Sigma A^T gives only three symmetric equations for four unknown A components, and requires the same material population with negligible births/deaths/visibility changes. Isotropic covariance is unchanged by arbitrary rigid rotation. Therefore statistics-only motion values are null; no full six-parameter recovery is claimed.

Official implementations and sources:

- [OpenPIV repository](https://github.com/OpenPIV/openpiv-python), [documentation](https://openpiv.readthedocs.io/en/latest/), [examples](https://github.com/OpenPIV/openpiv-python-examples). Installed **0.25.3, GPLv3**. Actual API has no Rust/backend selector; Rust comparison NOT_RUN for this version. We used `extended_search_area_piv`, `first_pass`, `multipass_img_deform`, validation and filters. Final replaced invalid vectors are explicitly excluded, not presented as measurements. The newer repository's speed claims are not our measurements.
- [pycpd](https://github.com/siavashk/pycpd), [API](https://siavashk.github.io/pycpd/), [Coherent Point Drift paper](https://arxiv.org/abs/0905.2635). Installed **2.0.0, MIT**. Official `AffineRegistration`, source previous frame Y, target current frame X, row convention Y@B+t; our column convention F=B.T. Three initial rotations diagnose ambiguity; convergence does not establish truth.
- [SEA-RAFT official code](https://github.com/princeton-vl/SEA-RAFT), [paper](https://arxiv.org/abs/2405.14793), [official small checkpoint](https://huggingface.co/MemorySlices/Tartan-C-T-TSKH-spring540x960-S). Code and weight model card BSD-3-Clause. Only this AI model was deployed.
- [WAFT official code](https://github.com/princeton-vl/WAFT) was reviewed as an alternative; not installed or run. SEA-RAFT-S was selected for bounded two-frame inference and lower setup burden on this 8 GB GPU, not because natural-image accuracy implies microscope accuracy.
- [SciPy least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html), [matrix logarithm](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.logm.html) supply the fit and continuous-time conversion.

Dependency license texts/pins are retained under `evidence/licenses/`; third-party code and large weights are not tracked in the research Git repository.
