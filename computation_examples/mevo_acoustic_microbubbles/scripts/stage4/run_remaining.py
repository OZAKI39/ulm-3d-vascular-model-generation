"""Finite ordered Stage 4 continuation; stop on error, never restart old jobs."""
import argparse
import json
import subprocess
import sys
import time
from scripts.stage4.common import ROOT,new_run,write,finish


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-run",required=True,help="Exact primary batch identity to wait for")
    parser.add_argument("--start",choices=["H","I","J","K","audit","report"],default="H")
    parser.add_argument("--secondary-certified-wall",action="store_true",help="Use verified exact wall-projection certificates in Stage I")
    args = parser.parse_args()
    if args.start in {"H","I"} and (ROOT/"results/stage4/performance_refactor/gold_baseline/GOLD_BASELINE_MANIFEST.json").exists():
        raise RuntimeError("Stage 4I gold is sealed. Do not rerun H/I; use the gated Stage 4I.5 completion workflow.")
    out = new_run("execution",vars(args))
    steps = [("H","experiments.py",["seeds"]),("H","statistics.py",[]),
             ("I","experiments.py",["secondary"]+(["--certified-wall"] if args.secondary_certified_wall else [])),("J","convergence.py",[]),
             ("K","history_branch.py",[]),("audit","final_audit.py",[]),("report","finalize.py",[])]
    state_path = ROOT/"results/stage4/primary_only/STATE.json"
    deadline = time.monotonic()+7200
    print("Waiting for exact Stage 4G batch: "+args.primary_run,flush=True)
    while True:
        if state_path.exists():
            state = json.loads(state_path.read_text())
            if state.get("run")==args.primary_run:
                if state["status"] not in ["PASS","PASS_WITH_SCOPE"]:
                    raise RuntimeError("Specified Stage 4G batch did not pass")
                break
        if time.monotonic()>deadline:
            raise TimeoutError("Finite wait for the specified Stage 4G completion expired")
        time.sleep(10)
    first = next(i for i,step in enumerate(steps) if step[0]==args.start)
    completed = []
    for stage,script,arguments in steps[first:]:
        command = [sys.executable,"-u",str(ROOT/"scripts/stage4"/script),*arguments]
        progress = dict(status="RUNNING",stage=stage,script=script,completed=completed,
                        primary_run=args.primary_run,execution=str(out.relative_to(ROOT)))
        write(out/"PROGRESS.json",progress)
        write(ROOT/"results/stage4/CURRENT_EXECUTION_STATE.json",progress)
        print("START_STAGE "+json.dumps(progress),flush=True)
        start = time.perf_counter()
        result = subprocess.run(command,cwd=ROOT)
        completed.append(dict(stage=stage,script=script,exit_code=result.returncode,elapsed_s=time.perf_counter()-start))
        if result.returncode:
            failure = dict(status="FAIL",completed=completed,failed_stage=stage,
                           scope="Finite ordered continuation stopped; partial scientific outputs remain unaccepted")
            write(ROOT/"results/stage4/CURRENT_EXECUTION_STATE.json",dict(failure,execution=str(out.relative_to(ROOT))))
            finish(out,failure)
            raise SystemExit(result.returncode)
    acceptance = json.loads((ROOT/"STAGE4_ACCEPTANCE.json").read_text())
    final = dict(status="COMPLETE",completed=completed,scientific_status=acceptance["stage_state"],execution=str(out.relative_to(ROOT)))
    write(ROOT/"results/stage4/CURRENT_EXECUTION_STATE.json",final)
    finish(out,final)


if __name__=="__main__":
    main()
