"""Confirm promising exact candidates without millions of component-timer calls."""
import json
import os
import platform
import subprocess
import sys

from scripts.stage4i5.common import BASE,ROOT,new_run,finish,write,sha
from scripts.stage4i5.compare import comparison


def completed_replay(backend,segment):
    """Reuse only closed evidence with identical replay sources and environment."""
    import importlib.metadata

    packages={name:importlib.metadata.version(name) for name in ["numpy","scipy","vtk","h5py"]}
    threads={name:os.environ.get(name) for name in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"]}
    matches=[]
    for path in sorted((BASE/"timing_confirmation").glob("*/SUMMARY.json")):
        row=json.loads(path.read_text())
        if (row.get("status")!="COMPLETE" or row.get("backend")!=backend
                or row.get("segment")!=segment or row.get("component_timers") is not False):
            continue
        identity=json.loads((path.parent/"IDENTITY.json").read_text())
        environment=identity["environment"]
        sources=identity["source_sha256"]
        relevant={name:digest for name,digest in sources.items()
                  if name.startswith("src/") or name in ["scripts/stage4i5/replay.py",
                      "scripts/stage4i5/variants.py","scripts/stage4/geometry_certificate.py"]}
        if (environment["python"]==platform.python_version() and environment["packages"]==packages
                and environment["threads"]==threads
                and identity["config"]["protocol_sha256"]==sha(BASE/"protocol/PROTOCOL.json")
                and all(sha(ROOT/name)==digest for name,digest in relevant.items())):
            matches.append(path.parent)
    return matches[-1] if matches else None


def confirm(candidates):
    """The original full gold is never rerun; only the four declared short windows."""
    protocol=json.loads((BASE/"protocol/PROTOCOL.json").read_text())
    out=new_run("timing_confirmation_report",dict(candidates=[r["backend"] for r in candidates],
        selection_margin=.05,scope="Pre-adoption timing control; same saved-state windows and physics"))
    runs={}
    for backend in ["P0",*[r["backend"] for r in candidates]]:
        for segment in protocol["segments"]:
            path=completed_replay(backend,segment)
            if path is None:
                write(out/"PROGRESS.json",dict(backend=backend,segment=segment["id"],status="RUNNING"))
                subprocess.run([sys.executable,"-u","-m","scripts.stage4i5.replay",
                    "--backend",backend,"--segments",segment["id"],"--no-component-timers"],
                    cwd=ROOT,check=True)
                path=completed_replay(backend,segment)
            if path is None:
                raise RuntimeError("No matching closed timing replay: "+backend+"/"+segment["id"])
            runs[(backend,segment["id"])]=path
    baseline=[]
    for evidence in candidates[0]["regression"]:
        row=comparison(ROOT/evidence["reference"],runs[("P0",evidence["segment"])])
        baseline.append(row)
    if not all(row["strict_P1_gate"]=="PASS" and row["bitwise_positions"] for row in baseline):
        write(out/"FAILURE.json",dict(reason="Removing timers changed P0 states",rows=baseline))
        raise RuntimeError("Minimal-timer P0 differs from instrumented P0")
    results=[]
    for candidate in candidates:
        backend=candidate["backend"]
        rows=[]
        timer_invariance=[]
        for evidence in candidate["regression"]:
            segment=evidence["segment"]
            rows.append(comparison(runs[("P0",segment)],runs[(backend,segment)]))
            timer_invariance.append(comparison(ROOT/evidence["candidate"],runs[(backend,segment)]))
        correct=all(row["strict_P1_gate"]=="PASS" for row in rows+timer_invariance)
        speed=sum(row["baseline_wall_s"] for row in rows)/sum(row["candidate_wall_s"] for row in rows)
        results.append(dict(backend=backend,strict_correctness=correct,aggregate_speedup=speed,
            regression=rows,timer_removal_regression=timer_invariance))
        write(out/"PARTIAL_RESULTS.json",results)
    return finish(out,dict(status="COMPLETE",candidates=results,baseline_timer_removal=baseline,
        scope="One serial measurement per unchanged window/backend; not a full 6 s optimized run. Component wrappers removed, common accepted-step diagnostics retained; no cProfile. Separate denominators from component-timed plots."))
