"""Shadow the compiled CPU contact projection on actual gold-derived states."""
import json
import time

import numpy as np
import pyvista as pv

from mevo_mb.bubbles.contact import project_contacts
from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep
from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4.geometry_certificate import install
from scripts.stage4i5.common import BASE,GOLD,ROOT,new_run,finish


def main():
    out=new_run("p1_vtk_contact",dict(candidate="P1_ordered",scope="Projection-only shadow; no trajectory acceptance"))
    compiled=OrderedContactSweep(ROOT/".optional/ordered_contact")
    rows=[]
    for segment in json.loads((BASE/"protocol/PROTOCOL.json").read_text())["segments"]:
        config=json.loads((GOLD/"cases"/segment["case"]/"RUN_IDENTITY.json").read_text())["config"]
        relative=config["geometry"]["mesh_surface"].split("/workspace/mevo_acoustic_microbubbles/")[-1]
        surface=pv.read(GOLD/"inputs"/relative)
        reference,candidate=[BifurcationGeometry(config["geometry"]) for _ in range(2)]
        for geometry in [reference,candidate]:
            geometry.bind_mesh_surface(surface)
            install(geometry,per_particle=True)
        with np.load(BASE/"protocol"/(segment["id"]+".npz")) as state:
            mask=state["active"]
            x,r,v=state["positions_m"][mask],state["radii_m"][mask],state["total_velocity_m_s"][mask]
        dt=min(config["particles"]["dt_s"],float(np.min(
            config["particles"]["max_displacement_per_radius"]*r/np.maximum(np.linalg.norm(v,axis=1),1e-300))))
        for fraction in [0.,.1,.25,.5,1.]:
            positions=x+fraction*dt*v
            args=(config["particles"]["contact_tolerance_m"],config["particles"]["contact_max_iterations"],"ckdtree")
            begin=time.perf_counter()
            a=project_contacts(positions,r,reference,*args)
            reference_s=time.perf_counter()-begin
            begin=time.perf_counter()
            b=compiled.project(positions,r,candidate,*args)
            candidate_s=time.perf_counter()-begin
            identical=bool(np.array_equal(a[0],b[0]) and np.array_equal(a[1],b[1]) and a[2:]==b[2:])
            rows.append(dict(segment=segment["id"],advance_fraction=fraction,advance_dt_s=fraction*dt,
                bitwise_and_events_identical=identical,max_position_error_m=float(np.max(np.linalg.norm(a[0]-b[0],axis=1))),
                iterations=a[3],contacts=len(a[2]),reference_s=reference_s,candidate_s=candidate_s))
    finish(out,dict(status="PASS" if all(r["bitwise_and_events_identical"] for r in rows) else "FAIL",
        candidate="P1_ordered",test="gold_projection_shadow",rows=rows,compiler_identity=compiled.identity,
        aggregate_projection_speedup=sum(r["reference_s"] for r in rows)/sum(r["candidate_s"] for r in rows),
        scope="Same 20 gold-derived inputs and original wall certificate, candidate selection, event order, tolerances and iteration bound. Full S1/S2 replay still required."))


if __name__=="__main__":
    main()
