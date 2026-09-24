"""Restore exact saved accepted states and the original binary interval schedule.

No historical state is recomputed. Rejected trials have no accepted physical
state and may be retried. Pending right-child endpoints/depths are reconstructed
from the final accepted binary leaf, retaining the original nominal dt grid.
"""
from collections import Counter
from dataclasses import replace
from types import MethodType
import numpy as np
from .particle6_stepper import bind_query_dependency
from .physical_time_refinement import refine_interval
from .lammps_state import BridgeParticle
from .particle_shapes import Sphere


def pending_intervals(samples,dt):
    last=float(samples[-1,0]);nearest=int(round(last/dt))
    if last==nearest*dt:return [],nearest+1
    nominal=int(np.floor(last/dt));begin=nominal*dt;end=(nominal+1)*dt
    leaf_begin=float(samples[-2,0]);leaf_end=last;depth=0;rights=[]
    while (begin,end)!=(leaf_begin,leaf_end):
        if depth>=48:raise ValueError('Saved accepted interval is not on original binary time tree')
        mid=begin+(end-begin)/2
        if leaf_end<=mid:
            rights.append((mid,end,depth+1));end=mid
        elif leaf_begin>=mid:begin=mid
        else:raise ValueError('Saved interval crosses an unsaved binary split')
        depth+=1
    pending=list(reversed(rights))
    if not pending or pending[0][0]!=last or pending[-1][1]!=(nominal+1)*dt:
        raise ValueError('Saved binary interval remainder not contiguous')
    if any(left[1]!=right[0] for left,right in zip(pending,pending[1:])):raise ValueError('Pending interval gap')
    return pending,nominal+2


def particle_at_saved_endpoint(meta,samples):
    row=samples[-1]
    return BridgeParticle.from_shape(meta['particle_id'],Sphere(row[1:4],meta['radius_m']),
        q=row[10:14],velocity=row[4:7],omega=row[7:10])


def restore_stepper_recording(stepper,meta,samples):
    stepper.samples=samples.tolist()  # float64 → Python float → float64 is exact
    stepper.accepted_count=meta['accepted_steps'];stepper.rejected_count=meta['rejected_trials']
    stepper.maximum_depth=meta['maximum_refinement_depth'];stepper.minimum_gap=meta['minimum_original_wall_gap_m']
    stepper.maximum_position_identity_error=meta['maximum_position_identity_error_m']
    stepper.states=Counter(meta['nearfield_states'])
    if not np.array_equal(stepper.read()[0].position,samples[-1,1:4]) or stepper.time_s!=samples[-1,0]:
        raise ValueError('Checkpoint stepper did not restore exact saved endpoint')


def advance_pending_interval(stepper,end,depth):
    # Keep the identical cached P6.5 trial and original P3 recursive arithmetic.
    # Only restore the depth of a previously pending right child.
    original=stepper.advance_cached
    def restored_refine(state,requested_end,trial,*,ledger):
        return refine_interval(state,requested_end,trial,depth=depth,ledger=ledger)
    function=bind_query_dependency(original.__func__,{'refine_interval':restored_refine})
    stepper.advance_cached=MethodType(function,stepper)
    try:return stepper.step_to(end)
    finally:stepper.advance_cached=original


def preceding_nominal_state(samples,dt):
    rows=[row for row in samples if row[0]==round(float(row[0])/dt)*dt]
    if not rows:raise ValueError('No original nominal state')
    unchanged=0
    for j in range(len(rows)-1,0,-1):
        if np.array_equal(rows[j][1:4],rows[j-1][1:4]):unchanged+=1
        else:break
    return np.array(rows[-1][1:4]),unchanged
