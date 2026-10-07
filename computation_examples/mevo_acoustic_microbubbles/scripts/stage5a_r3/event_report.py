"""Evidence-linked B2 root cause; no historical gate mutation."""

from scripts.stage5a_r3.common import BASE, read, write, sha


def main():
    evidence = read(BASE / "final/WALL_CERTIFICATE_AUDIT.json")
    assert all(
        r["cold_gap_difference_max_m"] == 0 and r["cold_normal_difference_max"] == 0
        for r in evidence["normal_switch"]
    )
    assert all(r["complete"] for r in evidence["rows"])
    trials = []
    for variant in ["nominal", "refined"]:
        source = BASE / "diagnosis/reproduction_v1" / variant / "RESULT.json"
        result = read(source)
        for r in result["trials"]:
            particle = next(p for p in r["exact_geometry"] if p["particle_id"] == 152)
            trials.append(
                dict(
                    variant=variant,
                    dt_s=r["dt_s"],
                    actual_time_change_s=r["actual_time_change_s"],
                    old_bound_m=r["old_safety"]["minimum_wall_swept_gap_bound_m"],
                    exact_gap_m=particle["minimum"]["gap_m"],
                    classification=particle["minimum"]["classification"],
                    particle=particle,
                    source=str(source),
                    source_sha256=sha(source),
                )
            )
    replays = []
    for path in sorted((BASE / "runs").glob("*local_v1_nominal/SUMMARY.json")):
        replays.append(read(path))
    for path in sorted((BASE / "runs").glob("*local_v1_refined/SUMMARY.json")):
        replays.append(read(path))
    hypotheses = [
        dict(
            id="H1",
            status="CONFIRMED",
            hypothesis="Global roundoff subtraction rejects a continuously safe terminal candidate",
            evidence="Exact rational finite-triangle squared distance passes both 1.52587890625 ps candidates; old lower bound fails",
        ),
        dict(
            id="H2",
            status="REFUTED",
            hypothesis="Removing the roundoff guard alone recovers B2",
            evidence="Actual full-system 1 ns and 50 ns candidates are exactly unsafe at the unchanged tolerance",
        ),
        dict(
            id="H3",
            status="CONFIRMED",
            hypothesis="Nearest-only WALL contact rows miss simultaneous finite-feature coverage",
            evidence="ID 152 contacts two finite face interiors (3178 and 13018). Adding their supporting contact rows permits both original full-system terminal states to advance 10 us with zero rejection/projection",
        ),
        dict(
            id="H4",
            status="REFUTED",
            hypothesis="The observed normal switch is an exact-coordinate cache artifact",
            evidence="Cold queries match cached gap/normals exactly. Independent exact nearest feature switches 3178 to 13018 on the same actual candidate",
        ),
        dict(
            id="H5",
            status="REFUTED",
            hypothesis="Terminal time/position cannot make a representable update",
            evidence="Original 1.52587890625 ps candidate changes absolute time by 1.5258905250448151 ps and moves ID 152 by about 3.37e-16 m. The next algorithmic request hits the 1 ps progress floor",
        ),
        dict(
            id="H6",
            status="UNRESOLVED",
            hypothesis="The complete 5 ms all-active motion has an independent error bound",
            evidence="Local convergence and cross-method comparisons do not establish a global trajectory error bound",
        ),
    ]
    report = dict(
        hypotheses=hypotheses,
        trials=trials,
        full_system_local_replays=replays,
        failure_semantics=dict(
            MINIMUM_PHYSICAL_PROGRESS="Algorithmic next dt below 1 ps after repeated rejected candidates; not proof that a binary64 update is impossible",
            NONLINEAR_NOT_CONVERGED="R2 S1 endpoint resistance alternates across actual nearest-feature normals; terminal diagnostics retained, no candidate promoted",
            SWEEP_NOT_CERTIFIED="Includes both safe-but-overconservative and genuinely unsafe candidates, distinguished by exact continuous geometry",
        ),
        tolerance_budget="Each contact row preserves its initial negative tolerance-band gap; the continuous certificate independently checks the original 1e-11 m bound. No repeated subtraction of tolerance, no clamping of raw gaps",
        unchanged=dict(
            radius_m=1.78e-6,
            active_B2=186,
            wall_mesh=True,
            nearest_wall_lubrication=True,
            activation_over_R=0.1,
            floor_over_R=0.001,
            projection_m=0,
        ),
        certificate_audit_sha256=sha(BASE / "final/WALL_CERTIFICATE_AUDIT.json"),
        repair_scope="Pure numerical finite-feature contact coverage and continuous predicate refinement; no added wall resistance or altered physical forcing",
    )
    write(BASE / "final/B2_EVENT_ROOT_CAUSE.json", report)


if __name__ == "__main__":
    main()
