"""Non-resetting budget ledger; distinguish measured work from conservative charges."""

import xml.etree.ElementTree as ET
from scripts.stage5a_r3.common import BASE, read, write, sha
from scripts.stage5a_r3.current import campaign


def main():
    rows = []
    # Leaf run summaries avoid counting aggregate LOCAL/EARLY/RESTART tables twice.
    for p in sorted((BASE / "runs").glob("*/SUMMARY.json")):
        if "formal_R3_" in p.parent.name:
            continue
        r = read(p)
        rows.append(
            dict(
                source=str(p),
                sha256=sha(p),
                category="development_trajectory",
                measured_s=r["wall_s"],
                allowance_s=50.0,
                reason="Body time plus conservative import/outer-wrapper allowance; allowance is not a timing measurement",
            )
        )
    for name, key in [
        ("diagnosis/COMMON_STATE_REPLAY_V1.json", "charged_experiment_wall_s"),
        ("diagnosis/LEGACY_ATTRIBUTION_V1.json", "charged_experiment_wall_s"),
        ("diagnosis/B2_LINEAGE_V1.json", "wall_s"),
        ("diagnosis/B3_LINEAGE_V1.json", "wall_s"),
        ("diagnosis/reproduction_v1/nominal/RESULT.json", "wall_s"),
        ("diagnosis/reproduction_v1/refined/RESULT.json", "wall_s"),
        ("final/WALL_CERTIFICATE_AUDIT.json", "wall_s"),
        ("benchmarks/LAYOUT_ISOLATED.json", "wall_s"),
        ("benchmarks/LAYOUT_ADOPTED_HORIZON_WITH_SYNC.json", "wall_s"),
        ("benchmarks/LAYOUT_ADOPTED_HORIZON.json", "wall_s"),
        ("benchmarks/LAYOUT_FULL_TRIAL.json", "wall_s"),
    ]:
        p = BASE / name
        if p.exists():
            rows.append(
                dict(
                    source=str(p),
                    sha256=sha(p),
                    category="development_diagnostic",
                    measured_s=read(p)[key],
                    allowance_s=50.0,
                    reason="Distinct diagnostic process/aggregate, not a trajectory-summary duplicate",
                )
            )
    p = BASE / "benchmarks/INDEPENDENT_THROUGHPUT.json"
    t = read(p)
    rows.append(
        dict(
            source=str(p),
            sha256=sha(p),
            category="independent_throughput",
            measured_s=None,
            charged_s=t["serial_wall_s"] + t["workers"] * t["parallel_wall_s"] + 50,
            reason="Parallel wall multiplied by worker count is a conservative aggregate charge, not single-task speedup",
        )
    )
    seen = set()
    for p in sorted((BASE / "environment").glob("JUNIT*.xml")):
        digest = sha(p)
        if digest in seen:
            continue
        seen.add(digest)
        root = ET.parse(p).getroot()
        suites = [root] if root.tag == "testsuite" else list(root.findall("testsuite"))
        rows.append(
            dict(
                source=str(p),
                sha256=digest,
                category="test_invocation",
                measured_s=sum(float(s.attrib.get("time", 0)) for s in suites),
                allowance_s=10.0,
                reason="Unique XML retained, including failed invocation; identical copied XML charged once",
            )
        )
    # Reserve for short exploratory checks without complete process latency logs.
    # This is deliberately included in budget use, never called measured runtime.
    rows.append(
        dict(
            category="unmetered_exploration_reserve",
            measured_s=None,
            charged_s=1800.0,
            reason="Conservative 30-minute budget reservation for earlier short probes/startups not individually metered; no performance claim is derived from it",
        )
    )
    for r in rows:
        if "charged_s" not in r:
            r["charged_s"] = r["measured_s"] + r["allowance_s"]
    total = sum(r["charged_s"] for r in rows)
    c = campaign()
    formal = dict(c["spent_s"])
    caps = dict(development=14400, new=14400, reference=28800, r2=5400, overall=63000)
    all_charged = total + sum(formal.values())
    checks = dict(
        development=total <= caps["development"],
        **{k: v <= caps[k] for k, v in formal.items()},
        overall=all_charged <= caps["overall"],
    )
    write(
        BASE / "final/BUDGET_LEDGER.json",
        dict(
            status="PASS" if all(checks.values()) else "FAIL",
            checks=checks,
            caps_s=caps,
            development_charged_s=total,
            development_rows=rows,
            formal_charged_s=formal,
            total_charged_s=all_charged,
            formal_prior_charges_carried_forward=c.get("prior_spent_s"),
            campaign_complete=c["complete"],
            limitations="Development body timing plus explicit conservative allowances is not a process-inclusive performance benchmark. V1 B2 refined lacks outer latency, so its run-body plus 50 s was carried forward as budget charge only. No budget is reset after revisions.",
        ),
    )
    print(dict(development_charged_s=total, formal_charged_s=formal, checks=checks))


if __name__ == "__main__":
    main()
