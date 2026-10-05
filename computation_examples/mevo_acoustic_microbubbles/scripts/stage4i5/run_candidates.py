"""Finite sequential I.5 candidate campaign, waiting for baseline timing to finish."""
import argparse
import json
import os
import subprocess
import sys
import time

from scripts.stage4i5.common import BASE,ROOT,new_run,finish,write


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--wait-for",required=True)
    args=parser.parse_args()
    out=new_run("candidate_execution",dict(wait_for=args.wait_for,concurrency=1))
    while True:
        status=subprocess.run(["supervisorctl","status",args.wait_for],text=True,capture_output=True).stdout
        write(out/"PROGRESS.json",dict(phase="WAITING_FOR_P0_BENCHMARK",supervisor_status=status.strip()))
        if "RUNNING" not in status and "STARTING" not in status:
            break
        time.sleep(15)
    complete=set()
    for path in (BASE/"p1_vtk_contact").glob("*/SUMMARY.json"):
        s=json.loads(path.read_text())
        if s.get("status")=="COMPLETE" and s.get("backend")=="P0" and not s.get("cprofile"):
            complete.add(s["segment"]["id"])
    if complete!={"B2_S1","B2_S2","B3_S1","B3_S2"}:
        raise RuntimeError("P0 short baseline incomplete; no candidate benchmark launched")
    operations=[]

    def run(module,arguments=(),optional=False,warp=False):
        env=os.environ.copy()
        if warp:
            env["PYTHONPATH"]=str(ROOT/".optional/warp-1.17.0")+":"+str(ROOT)
        command=[sys.executable,"-u","-m",module,*arguments]
        write(out/"PROGRESS.json",dict(phase=module,arguments=list(arguments),completed_operations=operations))
        code=subprocess.call(command,cwd=ROOT,env=env)
        operations.append(dict(command=command,exit_code=code,optional=optional))
        write(out/"OPERATIONS.json",operations)
        if code and not optional:
            raise RuntimeError("Required I.5 operation failed: "+module)
        return code

    run("scripts.stage4i5.profile_report")
    run("scripts.stage4i5.analyze_replays")
    run("scripts.stage4i5.wall_benchmark")
    run("scripts.stage4i5.neighbor_benchmark")
    for variant in ["P1_exact_cache","P1_neighbor","P1_combined","P1_static"]:
        run("scripts.stage4i5.replay",["--backend",variant,"--continue-on-failure"],optional=True)
        run("scripts.stage4i5.compare",["--backend",variant])
    if not run("scripts.stage4i5.install_warp",optional=True):
        run("scripts.stage4i5.wall_benchmark",["--backend","warp"],optional=True,warp=True)
        run("scripts.stage4i5.replay",["--backend","P2_warp","--continue-on-failure"],optional=True,warp=True)
        run("scripts.stage4i5.compare",["--backend","P2_warp"],optional=True)
    finish(out,dict(status="CANDIDATE_CAMPAIGN_FINISHED_REVIEW_REQUIRED",operations=operations,
        next="Review strict gates and limiter data; decide whether P3 is justified, then backend decision before J/K"))


if __name__=="__main__":
    main()
