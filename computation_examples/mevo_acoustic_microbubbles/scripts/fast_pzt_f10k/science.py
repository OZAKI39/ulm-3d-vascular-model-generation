"""Read-only scientific scope, exact original-wall admission, bounded saved states."""

import json
import subprocess
from time import perf_counter
import h5py
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import connected_components
from mevo_mb.fast_sim.vessel_geometry import VesselGeometry
from .common import ROOT, BASE, INPUT, BASELINE, identity, js, write, sha


class OriginalDomain:
    """Independent CPU triangle distance and original tetrahedral P1 interpolation."""

    def __init__(self):
        import pyvista as pv

        self.geometry = VesselGeometry(INPUT / "f4/labelled_lumen_m.vtp")
        self.grid = pv.read(INPUT / "result_121.vtu")
        assert set(self.grid.cells_dict) == {pv.CellType.TETRA}
        self.cells = np.asarray(self.grid.cells_dict[pv.CellType.TETRA])
        self.points = np.asarray(self.grid.points, float)
        self.values = np.asarray(self.grid["Velocity"], float)
        with h5py.File(INPUT / "f4/basis.h5") as h:
            # Actual archived schema, independent of a GPU field object.
            axes = [np.asarray(h["axes"][a + "_m"]) for a in ["x", "y", "z"]]
            self.lower = np.array([a[0] for a in axes])
            self.upper = np.array([a[-1] for a in axes])

    def sample(self, x, radius=1.78e-6):
        x = np.atleast_2d(x)
        distance, closest, normal, inside = self.geometry.query(x)
        tid = np.asarray(self.grid.find_containing_cell(x), int)
        velocity = np.full_like(x, np.nan)
        bary = np.full((len(x), 4), np.nan)
        valid = tid >= 0
        for i in np.flatnonzero(valid):
            tet = self.points[self.cells[tid[i]]]
            w = np.linalg.solve((tet[1:] - tet[0]).T, x[i] - tet[0])
            bary[i] = np.r_[1 - w.sum(), w]
            valid[i] = bool(bary[i].min() >= -1e-10 and bary[i].max() <= 1 + 1e-10)
            velocity[i] = bary[i] @ self.values[self.cells[tid[i]]]
        return dict(
            wall_gap_m=distance - radius,
            wall_closest_m=closest,
            outward_normal=normal,
            inside=inside,
            CFD_valid=valid,
            tetra_id=tid,
            barycentric=bary,
            CFD_velocity_m_s=velocity,
            field_valid=np.all((x >= self.lower + radius) & (x <= self.upper - radius), axis=1),
        )


def scope():
    reports = [
        "FAST_PZT_F5_FLOW_RELEASE_DEVELOPMENT.md",
        "FAST_PZT_F7_DEVELOPMENT.md",
        "FAST_PZT_F8_DEVELOPMENT.md",
        "FAST_PZT_F10A_DEVELOPMENT.md",
        "FAST_PZT_F10B_MODEL_SCOPE.md",
        "FAST_PZT_F10C_DEVELOPMENT.md",
        "FAST_PZT_F10D_DEVELOPMENT.md",
        "FAST_PZT_F10E_DEVELOPMENT.md",
        "FAST_PZT_F10F_DEVELOPMENT.md",
        "FAST_PZT_F10G_CONTACT_MODEL.md",
        "FAST_PZT_F10H_DEVELOPMENT.md",
        "FAST_PZT_F10I_DEVELOPMENT.md",
        "FAST_PZT_F10J_GPU_SOLVER_AUDIT.md",
        "FAST_PZT_F10J_PERFORMANCE.md",
    ]
    rows = []

    def row(model, task, status, evidence, limitation, files):
        rows.append(
            dict(
                model=model,
                task=task,
                status=status,
                evidence=evidence,
                limitation=limitation,
                sources={p: sha(ROOT / p) for p in files},
                training_ready=False,
            )
        )

    row(
        "M0",
        "Full-particle acoustics, CFD, contact and GPU batching",
        "NUMERICAL_CHECKED",
        "F10C B1/4/16/64 match original serial at fixed dt; F10D optimized postprocessing preserves metrics.",
        "Synchronous controls only. F10D/E distinct dense initial states reject. Engineering regression is not physical accuracy.",
        [
            "src/mevo_mb/fast_sim/f10c_engine.py",
            "src/mevo_mb/fast_sim/f10d_metrics.py",
            "docs/FAST_PZT_F10C_DEVELOPMENT.md",
        ],
    )
    row(
        "M0",
        "Interior 1kPa group response",
        "GROUP_MOTION_CHECKED",
        "F10A same-initial-state OFF/S/O10/5us gives resolved +Z paired centroid increment over10ms.",
        "One model and location; target navigation not established; F10K checks actual plane crossing.",
        ["docs/FAST_PZT_F10A_DEVELOPMENT.md", "results/fast_pzt/f10a/analysis/FINAL_GATE.json"],
    )
    row(
        "M0",
        "Internal morphology and graph topology",
        "MORPHOLOGY_LIMITED",
        "Complete saved IDs and accepted edge ledger; historical persistent changes in F5/F8 and preparations.",
        "Large dt-dependent individual positions, dense contact rejection; historical events are not verified B2 truth.",
        ["src/mevo_mb/fast_sim/morphodynamics.py", "docs/FAST_PZT_F10H_DEVELOPMENT.md"],
    )
    row(
        "B2",
        "Geometric unilateral contact, same discrete mobility",
        "NUMERICAL_CHECKED",
        "F10H independent two/four and real32; F10J corrected GPU true+recursive PCG and stable Armijo pass.",
        "Finite tested states, no universal convergence or unique physical equilibrium proof.",
        [
            "src/mevo_mb/fast_sim/f10h_contact.py",
            "src/mevo_mb/fast_sim/f10j_gpu.py",
            "results/fast_pzt/f10j/analysis/FINAL_GATE.json",
        ],
    )
    row(
        "B2",
        "Group motion and shape in one5ms state",
        "GROUP_MOTION_CHECKED",
        "F10H same-model10/5us: centroid0.048867nm, fixed-ID10.411515nm, Kabsch4.637528nm; shape thresholds pass.",
        "One O/1102 state, not positive navigation or arbitrary morphology validation.",
        ["docs/FAST_PZT_F10H_DEVELOPMENT.md", "results/fast_pzt/f10i/analysis/FINAL_GATE.json"],
    )
    row(
        "B2",
        "Strict edge timing",
        "MORPHOLOGY_LIMITED",
        "Original12/29 edge timings exceed10us; max130us. F10I tighter solve does not resolve timestep contribution.",
        "Keep original FAIL and order discrepancy; no threshold relaxation.",
        ["docs/FAST_PZT_F10I_EVENT_DIAGNOSTICS.md"],
    )
    row(
        "B2",
        "Persistent split and merge",
        "PERSISTENT_TOPOLOGY_UNTESTED",
        "All three original geometric B2 histories stay connected; no positive persistent event comparison beforeF10K.",
        "F10K candidate audit and conditional result are separate, do not retroactively validate historical events.",
        ["results/fast_pzt/f10j/analysis/TOPOLOGY_SCREEN.json"],
    )
    for model in ["M0", "B2"]:
        row(
            model,
            "Near-wall collective hydrodynamics",
            "NEAR_WALL_UNCALIBRATED",
            "Original triangle no-penetration is a geometric constraint.",
            "Not verified near-wall collective mobility; B2 interior adapter rejects any wall correction.",
            ["src/mevo_mb/fast_sim/f10h_diagnostic.py", "docs/FAST_PZT_F10G_CONTACT_MODEL.md"],
        )
        row(
            model,
            "Contact and near-field acoustics",
            "PHYSICS_UNCALIBRATED",
            "Original acoustic response and F2 exploratory closure retained.",
            "Local contact stiffness, lubrication cutoff and dense near-field force lack experimental calibration.",
            ["src/mevo_mb/fast_sim/nearfield.py", "docs/FAST_PZT_F10B_MODEL_SCOPE.md"],
        )
    write(
        BASE / "analysis/SCIENTIFIC_SCOPE.json",
        dict(
            rows=rows,
            reports={f"docs/{p}": sha(ROOT / "docs" / p) for p in reports},
            no_cross_model_truth=True,
            training_ready=False,
            identity=identity(),
        ),
    )


def inspect_state(path, frame, domain, event=None):
    r = 1.78e-6
    with h5py.File(path) as f:
        t = float(f["time_s"][frame])
        x = f["positions_m"][frame, 0]
        active = f["active" if "active" in f else "trajectory/active"][frame, 0]
        ids = f["particle_id"][:]
        radii = f["radii_m"][:]
        savedgap = f["wall_gap_m" if "wall_gap_m" in f else "trajectory/wall_gap_lower_bound_m"][frame, 0]
        cfg = (
            js(path.parent / "CONFIG.json")
            if (path.parent / "CONFIG.json").exists()
            else json.loads(f.attrs["config_json"])
        )
        controls = cfg["events"]
        applied = [e for e in controls if e["time_s"] <= t][-1]
        edge_t = f["edge_events/time_s"][:]
        edge_i = f["edge_events/batch_i_j_linked"][:]
        history_mask = edge_t <= t + 1e-12
        boundary_key = (
            "boundary_events_json" if "boundary_events_json" in f else "trajectory/boundary_events_json"
        )
        boundary = str(f[boundary_key][()]) if boundary_key in f else "NO_BOUNDARY_DATASET"
        source_identity = json.loads(f.attrs.get("identity_json", f.attrs.get("metadata_json", "{}")))
    spatial = domain.sample(x, r)
    d = squareform(pdist(x))
    np.fill_diagonal(d, np.inf)
    components, labels = connected_components(d <= 2.1 * r, directed=False)
    contacts = np.argwhere(np.triu(d < 2 * r, 1))
    nearest = float(d.min())
    reasons = []
    if not np.isfinite(x).all() or not active.all() or len(np.unique(ids)) != 32 or not np.all(radii == r):
        reasons.append("STATE_OR_ID_OR_RADIUS_INVALID")
    if not spatial["inside"].all() or spatial["wall_gap_m"].min() < -1e-12:
        reasons.append("FINITE_RADIUS_GEOMETRY_INVALID")
    if spatial["wall_gap_m"].min() < 10 * r:
        reasons.append("NEAR_WALL_OUTSIDE_INTERIOR_B2_SCOPE")
    if not spatial["CFD_valid"].all():
        reasons.append("CFD_DOMAIN_INVALID")
    if not spatial["field_valid"].all():
        reasons.append("STORED_ACOUSTIC_DOMAIN_INVALID")
    if max(0, 2 * r - nearest) > 17.8e-9:
        reasons.append("INITIAL_OVERLAP_BUDGET")
    name = path.parent.name
    dest = BASE / "analysis" / f"STATE_{name}_{frame}.npz"
    np.savez_compressed(
        dest,
        positions_m=x,
        particle_ids=ids,
        radii_m=radii,
        active=active,
        absolute_time_s=t,
        saved_gap_m=savedgap,
        labels=labels,
        active_contact_pairs=contacts,
        edge_history_time_s=edge_t[history_mask],
        edge_history_ids=edge_i[history_mask],
        **spatial,
    )
    return dict(
        source=str(path.relative_to(ROOT)),
        sha256=sha(path),
        frame=frame,
        absolute_start_s=t,
        state_file=str(dest.relative_to(ROOT)),
        state_sha256=sha(dest),
        particle_ids=ids,
        finite_positions=bool(np.isfinite(x).all()),
        all_active=bool(active.all()),
        radius_m=r,
        minimum_pair_distance_m=nearest,
        max_overlap_m=max(0, 2 * r - nearest),
        contacts=contacts,
        component_count=components,
        labels=labels,
        exact_wall_gap_min_m=float(spatial["wall_gap_m"].min()),
        wall_gap_recomputed_minus_saved_range_m=[
            float((spatial["wall_gap_m"] - savedgap).min()),
            float((spatial["wall_gap_m"] - savedgap).max()),
        ],
        exact_lumen_inside=bool(spatial["inside"].all()),
        CFD_valid=bool(spatial["CFD_valid"].all()),
        mean_CFD_velocity_m_s=spatial["CFD_velocity_m_s"].mean(0),
        field_valid=bool(spatial["field_valid"].all()),
        boundary_events=boundary,
        controls=controls,
        applied_control=applied,
        source_identity=source_identity,
        historical_event=event,
        edge_history_count=int(history_mask.sum()),
        restart=dict(
            positions_ids_radii=True,
            clock_absolute_offset=t,
            original_control_recoverable=True,
            physical_velocity_state="Overdamped Markov update; no inertial velocity, warm-start or tangential spring memory.",
            caches="Fresh coordinate-derived acoustic/CFD/wall caches; original field/flow update stride1; independently check M0 one-save-interval reconstruction before B2.",
            topology="M0 past ledger archived. New B2 graph initialized from actual positions; past dwell is left-censored and is not inherited as B2 history.",
            exact_B2_checkpoint=False,
        ),
        admission="B2_TOPOLOGY_CASE_NOT_ADMITTED"
        if reasons
        else "GEOMETRY_AND_STATE_ELIGIBLE_PENDING_GPU_RESTART_CHECK",
        reasons=reasons,
        training_ready=False,
    )


def main():
    started = perf_counter()
    (BASE / "analysis").mkdir(parents=True, exist_ok=True)
    old = ROOT / "results/fast_pzt/f10j"
    manifest = js(old / "ARTIFACT_MANIFEST.json")["files"]
    for name, v in manifest.items():
        assert sha(old / name) == v["sha256"], name
    inputs = js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json")
    for name, v in inputs.items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == v["sha256"], name
    write(
        BASE / "BASELINE_AUDIT.json",
        dict(
            F10J_complete_HEAD=BASELINE,
            verified_artifacts=len(manifest),
            verified_original_inputs=len(inputs),
            manifest_sha256=sha(old / "ARTIFACT_MANIFEST.json"),
            worktrees=subprocess.check_output(["git", "worktree", "list"], cwd=ROOT, text=True),
            local_AGENTS="No applicable AGENTS.md in project or ancestors",
            history_preserved=True,
        ),
    )
    scope()
    domain = OriginalDomain()
    historic = js(old / "analysis/TOPOLOGY_SCREEN.json")
    candidate = js(old / "analysis/TOPOLOGY_CANDIDATE.json")
    first = inspect_state(ROOT / candidate["source"], 835, domain, candidate["historical_event"])
    write(BASE / "analysis/F5_ADMISSION.json", first)
    alternatives = []
    # This finite order was committed in protocol before inspecting any candidates.
    chosen_paths = [
        "results/fast_pzt/f5/science_v1/main/trajectory.h5",
        "results/fast_pzt/f8/y_formal_v1/Y0/trajectory.h5",
        "results/fast_pzt/f10d/campaign_v2_diversity/prepare_1101/trajectory.h5",
        "results/fast_pzt/f10d/campaign_v2_diversity/prepare_1102/trajectory.h5",
    ]
    for name in chosen_paths:
        row = next(r for r in historic["rows"] if r["path"] == name)
        events = [
            e
            for e in row["topology"]["events"]
            if e["persistence_class"] == "persistent"
            and (e["kind"] == "merge" if "f5/" in name or "prepare" in name else e["kind"] == "split")
        ]
        event = (
            events[-1]
            if "prepare" in name
            else (next(e for e in events if e["time_s"] < 0.3) if "f5/" in name else events[0])
        )
        with h5py.File(ROOT / name) as f:
            times = f["time_s"][:]
            frame = max(0, int(np.searchsorted(times, event["time_s"] - 0.00025, side="right")) - 1)
        alternatives.append(inspect_state(ROOT / name, frame, domain, event))
    viable = [v for v in alternatives if not v["reasons"] and v["component_count"] > 1]
    selected = viable[0] if viable else None
    write(
        BASE / "analysis/ALTERNATIVE_SCREEN.json",
        dict(
            rows=alternatives,
            examined_states=len(alternatives),
            existing_inventory_sha256=sha(old / "analysis/TOPOLOGY_SCREEN.json"),
            inventory=historic["summary"],
            incomplete_ledger_exclusions=historic["reconstruction_status_counts"],
            existing_same_B2_histories=[r["path"] for r in historic["rows"] if r.get("is_same_geometric_B2")],
            selected=selected,
            no_random_search=True,
            no_new_coordinates=True,
            training_ready=False,
        ),
    )
    write(
        BASE / "analysis/SCIENCE_EXECUTION.json", dict(identity=identity(), wall_s=perf_counter() - started)
    )
    print("F5:", first["admission"], first["exact_wall_gap_min_m"], first["reasons"])
    print(
        "Alternative:",
        None
        if selected is None
        else (selected["source"], selected["frame"], selected["exact_wall_gap_min_m"]),
    )


if __name__ == "__main__":
    main()
