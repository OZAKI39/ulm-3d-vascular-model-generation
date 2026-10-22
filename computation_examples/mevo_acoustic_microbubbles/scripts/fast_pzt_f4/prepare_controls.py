"""Pre-trajectory phase alignment using the measured original PZT basis, no force edits."""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from mevo_mb.fast_sim.models import Basis, ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import BASE, ROOT, js, write


def main():
    p = F2Physics()
    b = Basis(BASE / "inputs/basis.h5", p.frequency_Hz)
    response = ResponseSection(BASE / "inputs/response.h5", p)
    amplitude = 291.237687721856
    rows = {}
    f3 = js(ROOT / "results/fast_pzt/f2/mature_v1/FIELD_SELECTION.json")["baseline"]
    for region in ["parent", "junction"]:
        x = np.load(BASE / "inputs" / (region + "_initial.npz"))["positions_m"]
        center = x.mean(0)
        h = np.array([RegularGridInterpolator(b.axes, z)(center[None])[0] for z in b.db.pressure_basis])
        phase = -np.angle(h)
        baseline = dict(time_s=0, amplitude_Pa=[amplitude] * 4, phase_rad=phase.tolist(), enabled=[True] * 4)
        switched = dict(baseline, phase_rad=(phase + np.array([np.pi / 6, 0, 0, 0])).tolist())
        info = {}
        for name, event in [("F3_reference", f3), ("baseline", baseline), ("switched", switched)]:
            field = b.combine(event["amplitude_Pa"], event["phase_rad"], event["enabled"])[0]
            v = field.interpolate_complex_pressure(x)
            response.query(v)
            info[name] = dict(
                event=event,
                pressure_min_Pa=float(abs(v).min()),
                pressure_RMS_Pa=float(np.sqrt(np.mean(abs(v) ** 2))),
                pressure_max_Pa=float(abs(v).max()),
            )
        rows[region] = dict(center_m=center, basis_per_source_at_center=h, **info)
    protocol = js(ROOT / "configs/fast_pzt/f4/protocol.json")
    protocol.update(
        aggregation_duration_s=0.6,
        control_selection="Before all acoustic trajectories: conjugate basis phase at initial cloud center, four amplitudes at the inherited 291.2376877 Pa campaign cap; switch Y+ phase by +30 degrees without amplitude/enable change",
        control_cap_scope="Previous F3 conservative campaign cap; not a measured hardware rating",
    )
    write(ROOT / "configs/fast_pzt/f4/protocol.json", protocol)
    # JSON separates complex values; raw complex basis samples remain in NPZ.
    raw = {n: rows[n].pop("basis_per_source_at_center") for n in rows}
    np.savez_compressed(BASE / "inputs/phase_alignment_basis.npz", **raw)
    write(
        ROOT / "configs/fast_pzt/f4/controls.json",
        dict(
            regions=rows,
            selection_before_acoustic_trajectories=True,
            reason="Original F3 zero-phase pressure is only 124/50 Pa at the two real near-wall sites; use actual basis phase alignment within inherited amplitude cap, with a bounded 600 ms formation window",
            no_model_change=True,
            no_post_trajectory_tuning=True,
        ),
    )
    print(
        {
            n: {k: v[k]["pressure_RMS_Pa"] for k in ["F3_reference", "baseline", "switched"]}
            for n, v in rows.items()
        }
    )


if __name__ == "__main__":
    main()
