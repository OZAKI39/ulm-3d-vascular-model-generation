"""Versioned spatial observation schema and explicit single-group legacy adapter."""
import copy,time
from xcommon import NAMES,empty
from state import make_state
def validate_channel(state):
 assert state['schema_version']=='M4-channel-3.0'
 assert state['training_ready'] is False and state['hardware_control_executed'] is False
 w=state['windows'];n=len(w['x']);assert all(len(v)==n for v in w.values())
 for i in range(n):
  if not w['valid'][i]:assert w['vx'][i] is None and w['vy'][i] is None
  if w['valid'][i]:assert w['fresh'][i] and not w['stale'][i] and w['safe'][i] and state['channel_mask_valid']
  if w['last_update_pts'][i] is not None:assert w['last_update_pts'][i]<=state['source_PTS']+1e-12
 for g in state['groups']:
  assert g['physical_id'] is None
  for name in NAMES:
   if not g['six']['valid'][name]:assert g['six']['values'][name] is None
 return True
def to_v2(state,group_key=None):
 """No implicit global fit: a caller must select one spatial observation group."""
 group=next((g for g in state['groups'] if g['observation_key']==group_key),None)
 est=copy.deepcopy(group['six']) if group else empty('EXPLICIT_LOCAL_GROUP_REQUIRED','F10X')
 est['observability_status']='VALID_GLOBAL_MOTION' if all(est['valid'].values()) else 'NOT_OBSERVABLE'
 dist={'group_center_xy':group['center_global_xy'] if group else None,'visible_count':group['visible_count'] if group else None,'field_validity':{'group_center_xy':group is not None},'scope':'EXPLICIT_SELECTED_SPATIAL_GROUP_NOT_WHOLE_CHANNEL'}
 out=make_state(state['source'],state['session_id'],state['frame_id'],state['source_PTS'],state['host_receive_timestamp'],state['ROI'],dist,est,{'channel_mask_valid':state['channel_mask_valid']},state['time_window'])
 out['schema_version']='M4-collective-2.0';out['model_provenance']={'name':state['backend'],'weight_sha256':'NOT_APPLICABLE_NO_TRAINED_MODEL','calibration_version':'F10X_PIXEL_ONLY_UNCALIBRATED'}
 out['stale_or_discontinuous']=state['stale_or_discontinuous'];out['selected_observation_key']=group_key
 return out
def to_v1(state,group_key=None):
 out=to_v2(state,group_key);out['schema_version']='M4-collective-1.0';out['motion']['observability_status']='INTERVAL_AFFINE_ESTIMATE' if all(out['motion']['valid'].values()) else 'NOT_OBSERVABLE';return out
def expire(state,now):
 age=(now-state['host_receive_timestamp'])*1000
 if age<=100:return state
 out=copy.deepcopy(state)
 out['stale_or_discontinuous']=True;out['age_ms']=age
 w=out['windows'];n=len(w['x']);w['valid']=[False]*n;w['fresh']=[False]*n;w['stale']=[True]*n;w['vx']=[None]*n;w['vy']=[None]*n
 out['fresh_valid_group_count']=0;out['valid_window_fraction_channel']=0.
 out.setdefault('reason_dictionary',{})[8]='STALE_OBSERVATION'
 w['reason_code']=[8 if safe else 1 for safe in w.get('safe',[True]*n)]
 if 'age_ms' in w:w['age_ms']=[v+age if v is not None else None for v in w['age_ms']]
 for g in out['groups']:g['six']=empty('STALE_OBSERVATION','F10X');g['visible_count']=None;g['center_global_xy']=None
 return out
