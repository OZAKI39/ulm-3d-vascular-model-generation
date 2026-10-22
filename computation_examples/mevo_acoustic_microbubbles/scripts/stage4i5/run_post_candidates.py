"""Finite I.5 completion -> revised J -> K -> audit, gated at every transition."""
import argparse
import json
import subprocess
import sys
import time

from scripts.stage4i5.common import BASE,ROOT,new_run,finish,write


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--wait-for",required=True)
    args=parser.parse_args()
    out=new_run("completion_execution",dict(wait_for=args.wait_for,
        order=["P1_ordered","backend_decision","J","K","final_audit","final_report"]))
    operations=[]
    while True:
        check=subprocess.run(["supervisorctl","status",args.wait_for],text=True,capture_output=True)
        if check.returncode not in [0,3]:
            raise RuntimeError("Cannot determine owned candidate service status: "+check.stderr)
        write(out/"PROGRESS.json",dict(phase="WAITING_FOR_CANDIDATE_CAMPAIGN",service=check.stdout.strip()))
        if "RUNNING" not in check.stdout and "STARTING" not in check.stdout:
            break
        time.sleep(15)
    candidate=json.loads((BASE/"candidate_execution/STATE.json").read_text())
    if candidate["status"]!="CANDIDATE_CAMPAIGN_FINISHED_REVIEW_REQUIRED":
        raise RuntimeError("Previous candidate campaign did not reach its natural completion boundary")

    def run(phase,module,arguments=(),optional=False):
        command=[sys.executable,"-u","-m",module,*arguments]
        progress=dict(stage=phase,status="RUNNING",command=command,completed_operations=operations,
            execution=str(out.relative_to(ROOT)),gold_rerun=False)
        write(out/"PROGRESS.json",progress)
        write(ROOT/"results/stage4/CURRENT_EXECUTION_STATE.json",progress)
        code=subprocess.call(command,cwd=ROOT)
        operations.append(dict(stage=phase,command=command,returncode=code,optional=optional))
        write(out/"OPERATIONS.json",operations)
        if code and not optional:
            write(ROOT/"results/stage4/CURRENT_EXECUTION_STATE.json",dict(progress,status="FAILED",returncode=code))
            raise RuntimeError("Required phase failed: "+phase)
        return code

    # The additional bottleneck-driven CPU candidate has its own independent
    # bitwise shadow, then exactly the original predeclared four replay windows.
    if run("I.5/projection_shadow","scripts.stage4i5.ordered_sweep_benchmark",optional=True)==0:
        shadows=[]
        for path in sorted((BASE/"p1_vtk_contact").glob("*/SUMMARY.json")):
            row=json.loads(path.read_text())
            if row.get("test")=="gold_projection_shadow":
                shadows.append(row)
        if shadows and shadows[-1]["status"]=="PASS":
            run("I.5/ordered_replay","scripts.stage4i5.replay",["--backend","P1_ordered","--continue-on-failure"],optional=True)
    run("I.5/ordered_comparison","scripts.stage4i5.compare",["--backend","P1_ordered"])
    # Explicit rejection/incomplete report is still required if optional Warp
    # installation failed. It cannot silently disappear from the decision table.
    run("I.5/warp_comparison","scripts.stage4i5.compare",["--backend","P2_warp"])
    run("I.5/backend_decision","scripts.stage4i5.backend_decision")
    for phase,module in [("J","convergence"),("K","history_branch"),("audit","final_audit"),("report","finalize")]:
        run(phase,"scripts.stage4."+module)
        if phase=="audit":
            audit=json.loads((ROOT/"results/stage4/audit/final/STATE.json").read_text())
            if audit["status"]!="PASS":
                raise RuntimeError("Final audit did not pass; no final report acceptance")
    result=finish(out,dict(status="COMPUTATION_AND_REPORT_COMPLETE_BACKUP_REVIEW_PENDING",operations=operations))
    write(ROOT/"results/stage4/CURRENT_EXECUTION_STATE.json",result)


if __name__=="__main__":
    main()
