"""Ideal scheduled events around the unchanged F10Q Sequential QP engine."""
import os
import time
import resource
import traceback
import subprocess
from dataclasses import replace
import numpy as np
import h5py
from mevo_mb.control.scheduled_pzt import Schedule
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from scripts.straight_tube_exp1.f10q.run import Campaign as ReferenceCampaign
from scripts.straight_tube_exp1.f10o.analyze import metrics, topology_record
from .common import ROOT, BASE, OLD, Q, CONFIG, context, policy, identity, read, write, sha


def shorten(schedule, duration):
    segments=[]
    for s in schedule.segments:
        if s.start_s >= duration: break
        segments.append(replace(s, duration_s=min(s.duration_s,duration-s.start_s)))
    return Schedule(tuple(segments))


class Campaign(ReferenceCampaign):
    def __init__(self,out):
        self.started=time.perf_counter();self.out=out;self.config=read(CONFIG);self.p=policy()
        self.ident=identity();self.records=[]
        c=context();self.x0=c["x0"];self.physics=c["physics"];self.response=c["response"];self.cavity=c["cavity"]
        self.bases={"coarse":c["basis"],"fine":CachedWaterBasis(OLD/"cache/fine.h5")}
        self.states=dict(np.load(BASE/"INITIAL_STATES.npz"))
        self.prior_s=sum(read(p)["wall_s"] for p in BASE.glob("*/CAMPAIGN.json"))
        self.route_budget=self.config["route_budget_s"][os.environ["F10S_MODE"]]
        self.prior_count=sum(len(read(p).get("runs",[])) for p in BASE.glob("*/CAMPAIGN.json"))
        self.schedules={k:Schedule.from_record(v) for k,v in read(BASE/"RUN_SCHEDULES.json").items()}
        for name,digest in read(Q/"PREREGISTRATION.json")["input_sha256"].items():
            assert sha(OLD/name)==digest,name
        for name in ["INITIAL_STATES.npz","RUN_SCHEDULES.json"]:
            assert sha(BASE/name)==read(BASE/"PREREGISTRATION.json")["input_sha256"][str((BASE/name).relative_to(ROOT))]

    def run(self,name,key,*,duration=None,dt=1e-6,cache="coarse",initial="original"):
        source_key=key
        if self.prior_s+time.perf_counter()-self.started>self.config["remote_budget_s"] or self.prior_count+len(self.records)>=self.config["max_trajectories"]:
            raise RuntimeError("Registered global calculation budget reached")
        schedule=self.schedules[key]
        if duration is not None:schedule=shorten(schedule,duration)
        folder=self.out/"runs"/name;folder.mkdir(parents=True,exist_ok=False)
        write(folder/"INPUT.json",dict(schedule=schedule.record(),source=key,dt_s=dt,duration_s=schedule.duration_s,
              cache=cache,initial=initial,identity=self.ident,preregistration_sha256=sha(BASE/"PREREGISTRATION.json"),
              schedule_manifest_sha256=sha(BASE/"RUN_SCHEDULES.json")))
        begin=time.perf_counter();self.x0=self.states[initial].copy();sim=self.engine(cache)
        x=[self.x0.copy()];times=[0.];steps=[];extra=[];pairs=[];fields=[];iterations=[];events=[];step_times=[]
        reasons=[];last=None;active=None
        integration=time.perf_counter()
        try:
            for j,(t,h,segment) in enumerate(schedule.steps(dt)):
                if time.perf_counter()-begin>self.config["per_run_budget_s"]:
                    reasons=["PER_RUN_WALL_BUDGET"];break
                if self.prior_s+time.perf_counter()-self.started>self.config["remote_budget_s"] or time.perf_counter()-self.started>self.route_budget:
                    reasons=["GLOBAL_REMOTE_WALL_BUDGET"];break
                if active != segment.segment_id:
                    before=sim.fields(sim.x) if active is not None else {"pressure_Pa":np.zeros(32,complex)}
                    prior=None if active is None else next(s.record() for s in schedule.segments if s.segment_id==active)
                    sim.control(segment.event());after=sim.fields(sim.x)
                    events.append(dict(ideal_command_time_s=segment.start_s,effective_time_s=t,segment_id=segment.segment_id,
                                       before_control=prior,after_control=segment.record(),positions_m=sim.x.copy(),
                                       before_pressure_Pa=before["pressure_Pa"],after_pressure_Pa=after["pressure_Pa"],
                                       timing="IDEAL_SCHEDULED_SEGMENTS",hardware_executed=False))
                    active=segment.segment_id
                last=sim.step(h);r=last["result"]
                step_times.append([t,h,active]);pairs.append(r["pair_force_N"].copy())
                extra.append([r["active_count"],r["local_mobility_residual_m"],r["local_force_residual_N"],
                              float(abs(last["fields"]["pressure_Pa"]).min()),float(abs(last["fields"]["pressure_Pa"]).max()),
                              float(np.linalg.norm(last["fields"]["total_force_N"],axis=1).max())])
                iterations.extend(dict(step=j,**rec) for rec in r["iterations"])
                if j%max(1,round(1e-4/dt))==0 or not last["accepted"]:
                    fields.append(dict(time_s=t,**last["fields"],contact_force_N=r["contact_force_N"]))
                steps.append([last["accepted"],r["mobility_equation_residual_m"],r["geometric_overlap_m"].max(initial=0),
                              r["max_single_iteration_m"],r["max_total_contact_correction_m"],last["wall_correction_m"],
                              last["full_step_s"],r["total_s"],last["resistance_s"],r["qp_iterations"],
                              np.linalg.norm(r["prediction_m"]-x[-1],axis=1).sum(),np.linalg.norm(r["contact_displacement_m"],axis=1).sum(),
                              r["primal_violation_m"],r["dual_negative_N"],r["force_equation_residual_N"],
                              r["complementarity_equivalent_gap_m"],r["complementarity_N_m"],r["actual_complementarity_equivalent_gap_m"]])
                if not last["accepted"]:reasons=last["reasons"];break
                x.append(sim.x.copy());times.append(t+h)
            if last is not None:write(folder/"LAST_STEP.json",last)
        except Exception as ex:
            reasons=["EXCEPTION: "+repr(ex)];(folder/"exception.txt").write_text(traceback.format_exc())
        integration_s=time.perf_counter()-integration
        x=np.asarray(x); rows=np.array(steps).reshape(-1,18)
        with h5py.File(folder/"trajectory.h5","w") as h5:
            for key,value in dict(positions_m=x,time_s=times,particle_id=np.arange(32),step_metrics=rows,
                                  pair_force_N=pairs,extra_metrics=extra,step_time_duration_segment=step_times).items():
                h5.create_dataset(key,data=value,compression="gzip")
            seq=h5.create_group("sequential_iterations")
            if iterations:
                for key in ["step","index","positions_m","linearization_m","pair_force_N","active_count","local_residual_m", "local_force_residual_N", "true_residual_m","true_force_residual_N","true_complementarity_m","overlap_m","working_set_iterations"]:
                    seq.create_dataset(key,data=np.array([v[key] for v in iterations]),compression="gzip")
            ff=h5.create_group("force_samples")
            if fields:
                for key in fields[0]:ff.create_dataset(key,data=np.array([f[key] for f in fields]))
            if reasons and last is not None:h5.create_dataset("failure_proposal_m",data=last["result"]["positions_m"])
            h5.attrs["flow"]="STRICT_ZERO";h5.attrs["execution_commit"]=self.ident["execution_commit"]
            h5.attrs["schedule_json"]=__import__('json').dumps(schedule.record())
            h5["step_metrics"].attrs["columns"]="accepted,residual_m,overlap_m,max_iteration_m,total_contact_m,wall_m,full_step_s,solver_s,resistance_s,iterations,physical_path_sum_m,contact_displacement_sum_m,primal_m,dual_N,force_residual_N,complementarity_m,complementarity_N_m,actual_complementarity_m"
            h5["extra_metrics"].attrs["columns"]="active_count,local_mobility_residual_m,local_force_residual_N,min_pressure_Pa,max_pressure_Pa,max_total_force_N"
        write(folder/"SWITCHES.json",events)
        result=dict(name=name,source=source_key,status="PASS" if not reasons and abs(times[-1]-schedule.duration_s)<1e-13 else "FAIL",
                    reasons=reasons,dt_s=dt,cache=cache,initial=initial,requested_s=schedule.duration_s,completed_s=times[-1],
                    accepted_steps=len(x)-1,attempted_steps=len(rows),displacement_m=x[-1].mean(0)-x[0].mean(0),
                    endpoint_m=x[-1].mean(0),integration_s=integration_s,wall_s=time.perf_counter()-begin,
                    mean_step_s=float(rows[:,6].mean()) if len(rows) else None,trajectory_sha256=sha(folder/"trajectory.h5"),
                    rss_peak_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                    mode=schedule.segments[0].mode,segment_count=len(schedule.segments),model="F10S_IDEAL_SCHEDULES_FROZEN_F10Q_SQP",
                    labels=self.config["labels"],Y_scope="HOMOGENEOUS_WATER_PROXY_SPECIFIC",training_ready=False,hardware_executed=False)
        for key,column in [("max_residual_m",1),("max_overlap_m",2),("max_contact_iteration_m",3),("max_contact_total_m",4),("max_wall_correction_m",5),("max_primal_m",12),("max_dual_N",13),("max_force_residual_N",14),("max_complementarity_m",15),("max_complementarity_N_m",16),("max_actual_complementarity_m",17)]:
            result[key]=float(rows[:,column].max(initial=0))
        write(folder/"RESULT.json",result);self.records.append(result);write(self.out/"RUNS.json",self.records)
        print(name,result["status"],"t_ms",times[-1]*1e3,"COM_um",result["endpoint_m"]*1e6,"wall_s",result["wall_s"],reasons,flush=True)
        return result


def eligibility(path):
    r=read(path/"RESULT.json")
    with h5py.File(path/"trajectory.h5") as h:t=h["time_s"][:];x=h["positions_m"][:]
    m=metrics(t,x);top,*_=topology_record(t,x)
    change=float(abs(m["rg_m"][-1]/m["rg_m"][0]-1))
    return dict(basic=r["status"]=="PASS",relative_Rg_change=change,connected=top["max_components"]==1,
                displacement_m=r["displacement_m"],eligible=bool(r["status"]=="PASS" and change<=.05 and top["max_components"]==1))


def main():
    mode=os.environ["F10S_MODE"]; label=os.environ["F10S_ATTEMPT"]
    out=BASE/label; out.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter(); c=None
    write(out/"RESOURCE_BEFORE.json",dict(identity=identity(),gpu=subprocess.check_output(["nvidia-smi","--query-gpu=name,memory.used,utilization.gpu","--format=csv,noheader"],text=True)))
    try:
        c=Campaign(out)
        if mode=="A":
            c.checks()
            for suffix,kw in [("DT1",{}),("DT05",dict(dt=.5e-6)),("FINE",dict(cache="fine"))]:
                if c.run("X_PAIR_3MS_"+suffix,"X_PAIR_3MS",**kw)["status"]!="PASS": break
        elif mode=="B":
            for key,sign in [("X_MINUS_ONLY",1),("X_PLUS_ONLY",-1)]:
                screen=c.run(key+"_SHORT",key,duration=.0005)
                eligible=screen["status"]=="PASS" and sign*screen["displacement_m"][0]>0
                write(out/(key+"_SELECTION.json"),dict(eligible=bool(eligible),screen=screen,rule="Complete 0.5ms with original basic gates and measured expected direction"))
                if eligible:
                    for suffix,kw in [("DT1",{}),("DT05",dict(dt=.5e-6)),("FINE",dict(cache="fine"))]:
                        if c.run(key+"_"+suffix,key,**kw)["status"]!="PASS": break
        else: raise ValueError(mode)
        status="COMPLETED"
    except Exception:
        status="EXCEPTION"; (out/"EXCEPTION.txt").write_text(traceback.format_exc()); raise
    finally:
        write(out/"CAMPAIGN.json",dict(status=locals().get("status","EXCEPTION"),mode=mode,wall_s=time.perf_counter()-start,
            runs=c.records if c else [],identity=identity(),peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            cpu_user_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime,cpu_system_s=resource.getrusage(resource.RUSAGE_SELF).ru_stime))

if __name__=="__main__":main()
