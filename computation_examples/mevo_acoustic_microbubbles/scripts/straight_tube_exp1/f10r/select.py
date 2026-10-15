"""Freeze trajectory controls using only completed field diagnostics."""
import time
from datetime import datetime, timezone
import numpy as np
from mevo_mb.control.scheduled_pzt import Schedule, axis_pair
from .common import BASE, CONFIG, context, read, write, identity, sha
from .diagnostics import primary


def main():
    begin=time.perf_counter()
    assert not (BASE/"RUN_SCHEDULES.json").exists()
    plans=read(BASE/"REGISTERED_SCHEDULES.json")
    plans.update(read(BASE/"FIELD_SELECTED_SCHEDULES.json"))
    proposal=read(BASE/"diagnostics/PROPOSALS.json")["selected_pair"]
    target=6*np.pi*.001*1.78e-6*2e-7/.002
    # One deterministic continuous-amplitude refinement within the registered source/pressure limits.
    # It is frozen before ANY new trajectory, and the original coarse proposal is preserved.
    scale=min(np.sqrt(target/proposal["mean_primary_N"][0]),1200/proposal["pressure_max_Pa"],1500/proposal["amplitude_Pa"])
    segment=axis_pair(x_amplitude=proposal["amplitude_Pa"]*scale,x_delta_deg=proposal["delta_deg"])
    c=context();pr,g,f,v=primary(c,segment.event(),c["x0"])
    assert abs(pr).max()<=1200*(1+1e-12)
    plans["X_PAIR"]=Schedule((segment,)).record()
    write(BASE/"RUN_SCHEDULES.json",plans)
    write(BASE/"CONTROL_REGISTRATION.json",dict(utc=datetime.now(timezone.utc).isoformat(),identity=identity(),
          fixed_target_protocol_sha256=sha(CONFIG),schedule_sha256=sha(BASE/"RUN_SCHEDULES.json"),
          original_pair_proposal=proposal,refined_pair_source_Pa=proposal["amplitude_Pa"]*scale,
          refined_pair_initial_mean_primary_N=f.mean(0),refined_pair_pressure_max_Pa=float(abs(pr).max()),
          selection_note="One deterministic source-amplitude scaling from the field-only table, within preregistered 1500 Pa-equivalent and 1200 Pa initial-particle limits. No trajectory used; target unchanged.",
          wall_s=time.perf_counter()-begin,particle_trajectories_run=0))
    print("Frozen X-pair source",segment.amplitude_Pa_equivalent,"force pN",f.mean(0)*1e12)


if __name__=="__main__":main()
