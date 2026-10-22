"""Read saved Stage 4 progress without opening active HDF5 files or starting jobs."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
KINDS = {"H": "statistics/paired_runs", "I": "primary_secondary",
         "J": "convergence", "K": "history_branch", "audit": "audit/final", "report": "report"}


def read(path):
    return json.loads(path.read_text()) if path.is_file() else None


def case_protocol(folder):
    identity=read(folder/"RUN_IDENTITY.json") or {}
    config=identity.get("config",{})
    particles=config.get("particles",{})
    return dict(secondary=identity.get("secondary_enabled"),
        nominal_dt_s=particles.get("dt_s"),displacement_per_radius=particles.get("max_displacement_per_radius"),
        sampling=config.get("sampling",{}).get("sampling_policy"),
        refined_field=config.get("effective_acoustic_field",{}).get("basis_override"),
        classification_tolerance_m=config.get("classification_tolerance_m"))


def performance_snapshot(base):
    """Only JSON metadata: incomplete optional candidates remain visible failures."""
    root=base/"performance_refactor"
    completed=[]
    pending=[]
    failures=[]
    now=datetime.now(timezone.utc).timestamp()
    for kind in ["profiles","p1_vtk_contact","p2_warp","timing_confirmation"]:
        for folder in sorted((root/kind).glob("*")):
            if not folder.is_dir():
                continue
            summary=read(folder/"SUMMARY.json")
            failure=read(folder/"EXECUTION_FAILURE.json")
            if failure:
                failures.append(dict(run=str(folder.relative_to(ROOT)),record=failure))
            elif summary and "segment" in summary and summary["status"]=="COMPLETE":
                completed.append(dict(run=str(folder.relative_to(ROOT)),backend=summary["backend"],
                    segment=summary["segment"]["id"],status=summary["status"],
                    profiled=summary.get("cprofile",False),component_timers=summary.get("component_timers",True),
                    accepted_steps=summary["accepted_steps"]))
            elif not summary and (folder/"PROGRESS.json").exists():
                path=folder/"PROGRESS.json"
                progress=read(path)
                pending.append(dict(run=str(folder.relative_to(ROOT)),evidence_kind=kind,**progress,
                    progress_age_s=round(max(0.,now-path.stat().st_mtime),1)))
    return dict(completed_short_replays=completed,unfinished_saved_replays=pending,
        candidate_failures=failures,backend_decision=read(root/"backend_decision/STATE.json"))


def snapshot():
    base = ROOT / "results/stage4"
    execution = read(base / "CURRENT_EXECUTION_STATE.json")
    result = dict(execution=execution, scope="Saved progress only; RUNNING is not a live process check.")
    if not execution or execution.get("status") != "RUNNING":
        return result
    if str(execution.get("stage","")).startswith("I.5"):
        result["performance_refactor"]=performance_snapshot(base)
        return result
    kind = KINDS.get(execution.get("stage"))
    if not kind:
        return result
    runs = sorted(p.parent for p in (base / kind).glob("*/IDENTITY.json"))
    if not runs:
        return result
    folder = runs[-1]
    # Do not present a previous invocation's batch as the currently queued work.
    if folder.name < Path(execution["execution"]).name:
        result["batch"] = "No batch created in this invocation yet"
        return result
    rows = read(folder / "ROWS.json") or []
    if not rows and kind=="convergence":
        for identity_path in sorted(folder.glob("*/RUN_IDENTITY.json")):
            case=identity_path.parent
            summary=read(case/"particles/SUMMARY.json")
            branch=read(case/"BRANCH_METRICS.json")
            if summary and summary.get("status")=="COMPLETE" and branch and not (case/"TERMINATION_POLICY_VIOLATION.json").exists():
                identity=read(identity_path)
                rows.append(dict(case=identity["control_case"],seed=identity["particle_seed"],
                    run=str(case.relative_to(ROOT)),protocol=case_protocol(case),
                    **{key:branch[key] for key in ["N_A","N_B","N_unassigned","N_wall_loss"]}))
    active = []
    failure_records = []
    failure_paths=list(folder.glob("*/EXECUTION_FAILURE.json"))+list(folder.glob("*/TERMINATION_POLICY_VIOLATION.json"))
    for path in sorted(failure_paths):
        try:
            record = read(path)
        except json.JSONDecodeError:
            record = dict(read_warning="Worker failure JSON was being written; retry this read-only command.")
        failure_records.append(dict(run=str(path.parent.relative_to(ROOT)),record=record))
    now = datetime.now(timezone.utc).timestamp()
    for path in sorted(folder.glob("*/PROGRESS.json")):
        if (path.parent / "BRANCH_METRICS.json").exists():
            continue
        try:
            progress = read(path) or {}
        except json.JSONDecodeError:
            progress = dict(read_warning="Worker progress JSON was being rewritten; retry this read-only command.")
        active.append(dict(run=str(path.parent.relative_to(ROOT)),
                           time_s=progress.get("time_s"),active_particles=progress.get("active_particles"),
                           accepted_steps=progress.get("accepted_steps"),protocol=case_protocol(path.parent),
                           elapsed_s=progress.get("elapsed_s"),
                           progress_age_s=round(max(0.,now-path.stat().st_mtime),1),
                           read_warning=progress.get("read_warning")))
    result.update(batch=str(folder.relative_to(ROOT)),completed_case_count=len(rows),
                  completed_cases=[{k:row[k] for k in ["case","seed","run","protocol","N_A","N_B","N_unassigned","N_wall_loss"] if k in row} for row in rows],
                  unfinished_saved_progress=active,failed_cases=read(folder/"FAILED_CASES.json") or [],
                  worker_failure_records=failure_records,
                  batch_summary=read(folder/"SUMMARY.json"))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",action="store_true",help="Print complete saved evidence pointers as JSON")
    args = parser.parse_args()
    state = snapshot()
    if args.json:
        print(json.dumps(state,indent=2,allow_nan=False))
        return
    current = state["execution"] or {}
    print("Recorded execution:",current.get("status","NOT_STARTED"),current.get("stage",""),current.get("script",""))
    if "batch" in state:
        print("Batch:",state["batch"])
        print("Completed cases:",state.get("completed_case_count",0))
    for row in state.get("unfinished_saved_progress",[]):
        print(Path(row["run"]).name,"time_s=",row["time_s"],"active=",row["active_particles"],
              "accepted_steps=",row["accepted_steps"],"elapsed_s=",row["elapsed_s"],"progress_age_s=",row["progress_age_s"])
        print("  Protocol:",json.dumps(row["protocol"]))
        if row["read_warning"]:
            print(row["read_warning"])
    for row in state.get("failed_cases",[]):
        print("FAILED:",json.dumps(row))
    for row in state.get("worker_failure_records",[]):
        print("WORKER_FAILURE:",json.dumps(row))
    if "performance_refactor" in state:
        perf=state["performance_refactor"]
        print("Completed short replays:",len(perf["completed_short_replays"]))
        for row in perf["unfinished_saved_replays"]:
            print(row.get("backend"),row.get("segment"),"time_s=",row.get("time_s"),
                "steps=",row.get("accepted_steps"),"evidence_kind=",row.get("evidence_kind"),
                "progress_age_s=",row["progress_age_s"])
        for row in perf["candidate_failures"]:
            print("CANDIDATE_FAILURE:",row["run"],row["record"].get("error"))
    print(state["scope"])


if __name__ == "__main__":
    main()
