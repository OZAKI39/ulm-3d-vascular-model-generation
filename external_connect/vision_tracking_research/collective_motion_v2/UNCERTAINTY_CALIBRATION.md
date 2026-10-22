# Independent synthetic calibration

Weights and validation rejection thresholds are frozen before final-test access. The independent calibration pool has 1,024 triples and96 class0 source-frame ×case groups. Each group score is its maximum standardized absolute residual, separately per parameter. The quantile uses ceil((G+1)*nominal) in sorted group maxima. This is an explicit finite-sample rank construction; no unverified MAPIE API is substituted. Source frames share a recording and cases reuse material, so exchangeability and independent-experiment coverage are not established. This procedure is an empirical conservative group calibration, not a real-domain guarantee.

| Model | Nominal | Actual six marginal coverages (%) | Six mean full interval widths |
|---|---:|---|---|
| Direct6 gray | 0.9 | 98.41, 96.67, 99.32, 96.67, 98.26, 98.26 | 88.199, 54.566, 1.636, 2.509, 0.704, 0.498 |
| Direct6 gray | 0.95 | 99.39, 98.94, 100.00, 99.77, 99.62, 99.09 | 106.560, 93.216, 1.968, 3.319, 0.899, 0.634 |
| Direct6 + heat | 0.9 | 98.33, 97.73, 95.68, 99.09, 98.94, 98.41 | 78.627, 62.612, 1.025, 2.885, 0.894, 0.640 |
| Direct6 + heat | 0.95 | 99.47, 98.86, 99.24, 99.92, 100.00, 99.92 | 90.942, 73.595, 2.209, 3.580, 1.231, 1.101 |
| GeoMotion | 0.9 | 98.48, 97.42, 89.92, 98.71, 95.15, 95.15 | 72.208, 47.717, 0.129, 2.655, 0.310, 0.156 |
| GeoMotion | 0.95 | 99.24, 98.86, 99.17, 99.55, 98.56, 99.55 | 80.555, 69.628, 1.826, 2.812, 0.739, 0.741 |

Coverage is often above nominal because group maxima create broad intervals. Broad intervals cannot rescue an inaccurate point estimate or the failed risk gate. GeoMotion's nominal90% divergence coverage is89.92%, while most95% marginals are above98%. Mean widths are provided to expose informativeness, and condition-specific coverage is in `analysis.json` and Figure5. No accepted-AI interval coverage exists because the operational accepted denominator is zero.

Source-frame bootstrap has only four source frames and one recording. Its intervals are descriptive, not independent physical-experiment confidence intervals. Real microscopy uncertainty remains UNCALIBRATED; no90/95% real guarantee is asserted.
