# F10L straight-tube experiment model

This is an independent scene based on the complete F10K code at
`850fd2df3225351f4fe9564f39cbe72a38662b30`.
The newer public branch `sync/mevo-acoustic-f10k-20261013` at
`15c42e0e422ff5d1521e5cc8c12ede2688b3de59` is a curated export of that source,
not a newer physical-model result or an ancestor of the complete local history.
F10L does not overwrite any Y-tube geometry, field, CFD or historical result.

## Geometry and meanings

The provisional inner diameter is 1.6 mm. Demo length is 20 mm and demo wall
thickness is 0.1 mm, giving a 1.8 mm outside diameter. Unknown experimental
thickness stays null in the lab parameter register. The finite lumen has
explicit inner/outer surfaces and two axial endpoints. Default endpoints are
`open_fail_stop`: a crossing raises bit 32768 and retains the last valid
saved frame, rather than introducing reflection or silently deleting an ID.
The optional closed endpoint includes the actual bubble radius in its constraint.

The default local X axis is axial. Three 9 mm circular acoustic faces are at
(0, +3.9, 0), (0, −3.9, 0), and (+13, 0, 0) mm, with normals (0, −1, 0),
(0, +1, 0), and (−1, 0, 0). The 3 mm provisional gap is measured from each
face to the outer tube envelope, or to the nearest endpoint plane. The two
side faces are 7.8 mm apart; each side face is 8.5 mm from the end face.
These are finite-disc distances, checked by constrained convex geometry,
not just distances between source centers.

Tube center/rotation, endpoint sign, source positions and normals are configured.
Two local side-Y sources and one axial source are enforced; no local Z± source
is admitted. The unknown injection hardware cannot be cleared for interference
without its dimensions. Rendered source faces have no invented housing thickness.

## No-flow full-particle dynamics

`ZeroFlow` and `GPUZeroFlow` return exact zero velocity and gradient without a
file loader. The new run never constructs the old CFD or triangular Y-wall
provider. An open-file audit rejects historical Y data during the trajectory.
The analytic side constraint is radial distance ≤ lumen radius minus bubble
radius. Small corrections are recorded; a proposal exceeding the original
0.1R wall budget is rejected. Original M0 residual contact tolerance 0.01R and
single-pass contact limit 0.1R remain unchanged.

`StraightTubeSimulator` subclasses the existing full-particle batch/F6 path.
F2 exploratory pair feedback, primary and secondary Bjerknes forces, overdamped
Stokes velocity, local lubrication, ordered reductions and PBD remain active.
There is no rigid cluster, adhesion, coalescence, new soft-contact law or B2
production substitution. A narrow explicit F3 opt-in allows response-matched
physical parameters while retaining the frozen numerical/contact choices;
the original default rejects non-frozen physics as before.

The demo starts after injection, with 32 fixed-ID equal-radius 1.78 µm bubbles
in a seed-20261013 cloud of radius 40 µm centered at (50,50,50) µm.
No continuous injection is present. Other counts, cloud centers, spreads, seeds
and durations are configurable. The current M0 adapter explicitly rejects a
polydisperse request rather than silently replacing it.

OFF and THREE_PZT_ON are registered for 1 ms at 1 µs, with 10 µs saved frames.
ON has three 300 Pa peak *numerical source* inputs at zero phase.
Saved-frame pressure and force observables are recomputed at actual stored
positions with the same GPU kernels; they are not claimed to be every-substep
force histories. All trajectories retain physical times and particle IDs.

## Acoustic model and response limits

Exactly three new k-Wave disc basis fields are computed, never the old four
6 mm fields. The homogeneous-water demo uses c=1500 m/s, density 1000 kg/m³,
and f=2.5 MHz. Unknown tube material is not an invented impedance layer:
the tube constrains bubble centers mechanically, but the propagation medium
is uniform. Transparent additive pressure sources do not model a reflective
ceramic body, backing, electric drive or measured electromechanical transfer.

The grid is derived from the full tube/source envelope and a 34-cell BLI margin.
At h=0.1 mm the wavelength has six cells; this is not grid convergence.
The outside PML has 16 cells. Actual BLI occupancy must remain clear of the
physical-grid edge. The requested CFL is 0.25, with 24 samples/period, four
smooth ramp periods, two diagonal transit times and four recording periods.
The two recording halves must agree within relative L2 0.01.

The convention is Re[P exp(+iωt)]. The database stores complex pressure per
source Pa and the F6 derivative of the same trilinear complex pressure.
At grid nodes this is the positive-cell derivative (last cell at the maximum);
the derivative is generally discontinuous across cells. Runtime re-derives
the gradient from pressure instead of interpolating an unrelated nodal gradient.
This consistency does not improve physical grid resolution.

The archived APECSS table SHA-256 is
`5d53cb2ef42fa64591f45e3e7de5db106d6cbe74794c8a9abfa88a3354bec983`.
Its actual bounds are R=1.779–1.781 µm, f=2.499–2.501 MHz and
nonzero pressure **0.1 Pa–40000 Pa**. The task's shorthand “0.1–40 kPa”
must not be used to relabel the first knot: the HDF5 data and its SI units are
authoritative. Exact zero drive has the exact quiescent response. All other
out-of-domain radii, frequencies and pressures are rejected, not clipped.

Source peak Pa, computed local peak |P|, local RMS |P|/√2, voltage and power
are separate quantities. No source-voltage conversion has been calibrated.
Three independent complex controls do not establish arbitrary 3D navigation:
the symmetric geometry limits spatial authority, and acoustic wavelengths
also limit microscopic shape addressing.

## Explicit omissions and future interface

Not modeled/calibrated: tube material reflection/absorption and elastic shear,
PZT body/backing/electrical transfer, measured bubble shell/radius distribution,
collective nonlinear radial dynamics, calibrated near-field pair forces,
wall hydrodynamics, full collective Stokes contact, inertia/history/added mass,
streaming, buoyancy, residual injection flow, Brownian motion and two-way fluid
coupling. “No active flow” in the lab does not imply these motions are zero.

The lab register includes future inside/wall/outside material data and acoustic
boundaries. The scalar adapter explicitly remains homogeneous. Heterogeneous
maps and a suitable backend must be added from measured data; a fluid-equivalent
solid wall would still omit elastic shear.

Target point/region, tolerance, target time, shape descriptor, three amplitudes,
three phases, event times and success/failure are reserved interfaces only.
This round records local individual-source gradients/forces at one location.
It does not optimize control or claim arrival, sustained aggregation, splitting,
or a prescribed shape.

Persistent status: EXPERIMENTAL_GEOMETRY_PENDING,
MATERIAL_ACOUSTICS_UNCALIBRATED (lab interface alias:
MATERIAL_ACOUSTICS_NOT_CALIBRATED), SOURCE_PRESSURE_UNCALIBRATED,
PHYSICS_UNCALIBRATED; training_ready=false, experiment_matched=false.

## References inspected

- [k-Wave array documentation](https://k-wave.org/documentation/kWaveArray.php)
  was unavailable (HTTP 502) during this audit; the actual installed 0.6.2
  `kWaveArray.add_disc_element` implementation and repository adapter were read.
- [Unified API](https://k-wave-python.readthedocs.io/en/latest/get_started/new_api.html):
  Python/CuPy and C++ backends, external PML, scalar or mapped medium inputs.
- [Official heterogeneous-medium example](https://github.com/waltsims/k-wave-python/blob/master/examples/ivp_heterogeneous_medium.py):
  reference for future measured material maps, not implemented physics here.
- [APECSS](https://github.com/polycfd/apecss):
  the existing response table is reused only in its registered domain.
- [NVIDIA Warp](https://nvidia.github.io/warp/):
  the existing float64 full-particle GPU execution is retained.
