"""Versioned observation-only M4 adapter. Existing 1.0 implementation unchanged."""
import copy,time
from legacy import make_state,ConsumerGate,invalidate

def new_state(*args,model_metadata,**kwargs):
    state=make_state(*args,**kwargs);state['schema_version']='M4-collective-2.0';state['model_provenance']=model_metadata
    state['scientific_accuracy_verified']=False;return state

def to_v1(state):
    result=copy.deepcopy(state);result['schema_version']='M4-collective-1.0'
    status=result['motion']['observability_status'];valid=result['motion']['valid']
    if status=='NON_AFFINE_MOTION':
        for k in valid:valid[k]=False;result['motion']['values'][k]=None;result['motion']['reasons'][k]='NON_AFFINE_MOTION_V2'
    result['motion']['observability_status']='INTERVAL_AFFINE_ESTIMATE' if all(valid.values()) else ('PARTIAL' if any(valid.values()) else 'NOT_OBSERVABLE')
    result['original_observability_status_v2']=status
    return result

class VersionedGate:
    def __init__(self):self.gate=ConsumerGate();self.last=None;self.session=None;self.last_pts=None
    def accept(self,state,now=None):
        gap=self.session==state['session_id'] and self.last is not None and (state['frame_id']!=self.last+1 or (self.last_pts is not None and state['source_PTS']-self.last_pts>.15))
        original=copy.deepcopy(state)
        if gap:original=invalidate(original,'FRAME_GAP_OR_PTS_DISCONTINUITY')
        # The old consumer receives only old enum values via the explicit adapter.
        checked,row=self.gate.accept(to_v1(original),now)
        if not row['fresh'] or gap:original=invalidate(original,'FRAME_GAP_OR_PTS_DISCONTINUITY' if gap else row['reason'])
        self.last=state['frame_id'];self.last_pts=state['source_PTS'];self.session=state['session_id'];self.latest=original
        return original,{**row,'fresh':row['fresh'] and not gap}
    def current(self,now=None):
        v1=self.gate.current(now)
        if v1 is None:return None
        if v1['stale_or_discontinuous']:return invalidate(self.latest,'NO_NEW_FRAME_EXPIRED')
        return self.latest
