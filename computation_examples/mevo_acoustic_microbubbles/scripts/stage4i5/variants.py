"""Explicit opt-in candidates; no production default changes before validation."""
import json
import pyvista as pv

import mevo_mb.bubbles.contact as contact
import mevo_mb.bubbles.tracker as tracking
from mevo_mb.bubbles.neighbor_cache import ContactNeighborCache
from mevo_mb.coupling.wall_backends import ExactCoordinateWallCache, PersistentVTKWallBackend

ORIGINAL_CANDIDATES = contact.contact_candidates
ORIGINAL_PROJECTION = tracking.project_contacts


def validate_runtime_identity(runtime,decision):
    """Formal compiled execution must use the binary/BLAS that passed regression."""
    if runtime["backend"]!="P1_ordered":
        return
    from scripts.stage4i5.common import ROOT
    selected=next(row for row in decision["candidates"] if row["backend"]==runtime["backend"])
    if not selected["strict_correctness"] or len(selected["regression"])!=4:
        raise ValueError("Ordered contact backend lacks four passing regression segments")
    for row in selected["regression"]:
        proof=json.loads((ROOT/row["candidate"]/"SUMMARY.json").read_text())
        if runtime["ordered_contact_sweep"]!=proof["backend_info"]["ordered_contact_sweep"]:
            raise ValueError("Compiled contact runtime identity differs from the accepted replay binary/BLAS")


def install_variant(tracker, backend):
    options = {
        "P1_exact_cache": (True,False,False),
        "P1_neighbor": (False,True,False),
        "P1_combined": (True,True,False),
        "P1_static": (True,False,True),
        "P2_warp": (True,False,False),
        "P1_ordered": (False,False,False),
    }
    if backend not in options:
        raise ValueError("Unknown candidate backend: "+backend)
    cache,neighbor,static=options[backend]
    info=dict(backend=backend,physics_changed=False)
    cleanup=[]
    if backend=="P1_ordered":
        from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep
        from scripts.stage4i5.common import ROOT
        sweep=OrderedContactSweep(ROOT/".optional/ordered_contact")
        info["ordered_contact_sweep"]=sweep.identity
        tracking.project_contacts=sweep.project
        cleanup.append(lambda:setattr(tracking,"project_contacts",ORIGINAL_PROJECTION))
    if cache:
        wall=PersistentVTKWallBackend(pv.read(tracker.config["geometry"]["mesh_surface"])) if static else None
        if backend=="P2_warp":
            from mevo_mb.coupling.warp_wall import WarpWallBackend
            wall=WarpWallBackend(pv.read(tracker.config["geometry"]["mesh_surface"]))
        exact=ExactCoordinateWallCache(tracker.cylinder,query_backend=wall)
        tracker.cylinder.level_and_gradient=exact.level_and_gradient
        tracker._i5_wall_cache=exact
        tracker._i5_wall_backend=wall
    if neighbor:
        candidates=ContactNeighborCache(4*tracker.p["size_distribution"]["radius_m"],tracker.p["contact_tolerance_m"])
        tracker._i5_neighbor_cache=candidates

        def cached(x,radii,gap_m=0.,backend="brute",timings=None):
            if backend!="ckdtree":
                return ORIGINAL_CANDIDATES(x,radii,gap_m,backend,timings)
            return candidates.candidates(x,radii,gap_m,timings)

        contact.contact_candidates=cached
        cleanup.append(lambda:setattr(contact,"contact_candidates",ORIGINAL_CANDIDATES))
    if cleanup:
        def restore():
            for action in reversed(cleanup):
                action()
        tracker._i5_cleanup=restore
    return info
