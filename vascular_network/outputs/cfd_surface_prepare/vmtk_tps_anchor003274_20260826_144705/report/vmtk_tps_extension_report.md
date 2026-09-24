# VMTK TPS flow-extension failure report

- Final status: `VMTK_EXTENSION_GEOMETRY_FAILED`
- NEXT: `REVIEW VMTK TPS FAILURE`
- Formal candidate count: `1`
- Second run or parameter sweep: `false`
- VMTK runtime: `1.5.0` (`v1.5.0`, commit `30d5d7cb8e607d153c208a9d7d39c9feb7985476`)
- VTK runtime: `9.2.6`
- Official filter: `vtkvmtkPolyDataFlowExtensionsFilter`
- Interpolation: `thinplatespline`
- Extension mode: `centerlinedirection`
- Transition ratio: `0.5`
- Preserve cross-section shape: `false`
- Custom TPS: `false`
- Automatic fallback: `false`

## Mandatory geometry evidence

| Port | Length (um) | Length error | Area error | Direction dot | Interface P95/P99 (deg) |
|---|---:|---:|---:|---:|---:|
| cut_000 | 15.563496 | 1.9555% | 0.2527% | 0.993813 | 12.388 / 15.057 |
| cut_001 | 11.760789 | 0.2165% | 0.6596% | 0.997909 | 8.524 / 10.198 |
| cut_002 | 6.415595 | 44.3504% | 36.0879% | 0.549542 | 58.110 / 67.060 |
| terminal_000 | 12.655853 | 0.8434% | 0.1361% | 0.998068 | 11.933 / 12.682 |

The fixed direction threshold is `0.999`; it was not relaxed. The RAW surface has one component, four open profiles, zero nonmanifold edges, zero self-intersections, zero degenerate triangles, and zero extension collisions. Nevertheless, the mandatory direction test fails at all four ports, and `cut_002` also fails length and distal-area fidelity.

The VMTK interface P95/P99 is worse than the old custom reference at every port. Therefore the answer to whether this candidate clearly removes the old visible fold is `NO`. The intermediate remeshed/capped VTP is retained only as failure evidence and is not promoted to a CFD STL/VTP. Radius fidelity, core fidelity, pressure correction, final tagging, and review figures were not run after the mandatory RAW geometry failure.

Original Ultraliser STL/VTP and old-custom comparison STL/VTP hashes were rechecked and are unchanged.
