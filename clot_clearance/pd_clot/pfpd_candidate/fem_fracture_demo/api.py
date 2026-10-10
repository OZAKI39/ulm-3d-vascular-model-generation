"""Unit-explicit provider contract, independent of the downstream clot model."""
from pathlib import Path
import json,math
PROVIDERS={'USER_SPECIFIED','LITERATURE_GUIDED','LITERATURE_GUIDED_USER_CONFIG'}
FUTURE_PROVIDERS={'BUBBLE_REDUCED','CFD_REDUCED'}
DYNAMIC={'frequency','frequency_hz','omega','omega_rad_s','impact','fast_impact','bubble_collapse','pulse_duration_s','time_step_s','velocity_m_s','strain_rate_s_1'}
class RequestError(ValueError):
 def __init__(self,code,message):self.code=code;super().__init__(f'{code}: {message}')
def fail(code,message):raise RequestError(code,message)
def validate(r,output_dir=None):
 def walk(x):
  if isinstance(x,dict):
   if any(k.lower() in DYNAMIC or 'frequency' in k.lower() or 'collapse' in k.lower() for k in x):fail('DYNAMIC_MODEL_REQUIRED','A dynamic loading key requires an inertial model')
   for v in x.values():walk(v)
  elif isinstance(x,(list,tuple)):
   for v in x:walk(v)
  elif isinstance(x,(float,int)) and not isinstance(x,bool) and not math.isfinite(x):fail('INVALID_REQUEST','All numeric values must be finite')
  elif isinstance(x,str) and x.lower().replace('_',' ') in {'fast impact','bubble collapse'}:fail('DYNAMIC_MODEL_REQUIRED','Fast impact or bubble collapse requires an inertial model')
 walk(r)
 if r.get('loading_regime','QUASISTATIC')!='QUASISTATIC':fail('DYNAMIC_MODEL_REQUIRED','Only quasi-static loading is implemented')
 if r.get('provider_id') not in PROVIDERS:fail('PROVIDER_NOT_IMPLEMENTED' if r.get('provider_id') in FUTURE_PROVIDERS else 'UNKNOWN_PROVIDER','Provider identity must be explicit')
 for k in ['region','direction','source','units']:
  if not r.get(k):fail('INVALID_REQUEST',f'Missing {k}')
 if not isinstance(r['units'],dict):fail('INVALID_REQUEST','Units must be a mapping')
 d=r['direction']
 if not isinstance(d,list) or len(d)!=3 or not all(isinstance(v,(float,int)) for v in d) or sum(v*v for v in d)<=0:fail('INVALID_REQUEST','A nonzero 3-vector direction is required')
 c=r.get('load_control')
 if c=='DISPLACEMENT_CONTROL':
  for k in ['loading_path','max_steps']:
   if k not in r:fail('INVALID_REQUEST',f'Missing {k}')
  keys=[k for k in ['max_displacement_m','stretch_ratio'] if k in r]
  if len(keys)!=1:fail('INVALID_REQUEST','Declare one displacement limit or stretch ratio')
  key=keys[0];unit='m' if key=='max_displacement_m' else '1'
  if r['units'].get(key)!=unit or r['units'].get('loading_path')!=unit:fail('INVALID_UNITS',f'{key} and path require {unit}')
  path=r['loading_path'];limit=r[key]
  if not isinstance(path,list) or len(path)<2 or not isinstance(r['max_steps'],int) or r['max_steps']<len(path)-1:fail('INVALID_REQUEST','Loading path and positive step budget required')
  if any(b<=a for a,b in zip(path,path[1:])) or path[-1]>limit or path[0]!=(0. if key=='max_displacement_m' else 1.):fail('INVALID_REQUEST','Path must increase from the reference value within its declared limit')
 elif c=='TRACTION_CONTROL':
  for k in ['center_m','width_m','load_steps','max_amplitude','p_Pa','tau_Pa']:
   if k not in r:fail('INVALID_REQUEST',f'Missing {k}')
  for key,unit in {'center_m':'m','width_m':'m','p_Pa':'Pa','tau_Pa':'Pa','max_amplitude':'1','load_steps':'1'}.items():
   if r['units'].get(key)!=unit:fail('INVALID_UNITS',f'{key} requires {unit}')
  if len(r['center_m'])!=3 or r['width_m']<=0 or not isinstance(r['load_steps'],list) or not r['load_steps'] or max(r['load_steps'])>r['max_amplitude']:fail('INVALID_REQUEST','Invalid traction support or amplitude path')
 else:fail('INVALID_REQUEST','Unknown control mode')
 if output_dir is not None and Path(output_dir).exists():fail('OUTPUT_EXISTS','Existing case directories are immutable, including successful cases')
 return r
def load_provider(r):
 """Both control modes normalize to a fixed, provider-neutral boundary payload."""
 validate(r);keys=['provider_id','load_control','region','direction','source','units']
 payload={k:r[k] for k in keys};payload['model_time']='QUASISTATIC_PATH_PARAMETER_NOT_PHYSICAL_TIME'
 if r['load_control']=='DISPLACEMENT_CONTROL':payload['prescribed_kinematics']={k:r[k] for k in ['max_displacement_m','stretch_ratio','loading_path','max_steps'] if k in r};payload['z_load']=None
 else:payload['z_load']={k:r[k] for k in ['center_m','width_m','p_Pa','tau_Pa','load_steps','max_amplitude']};payload['prescribed_kinematics']=None
 return payload
def create_case(r,output_dir):
 validate(r,output_dir);p=Path(output_dir);p.mkdir(parents=True,exist_ok=False);(p/'REQUEST.json').write_text(json.dumps(r,indent=2,allow_nan=False));(p/'LOAD_PROVIDER.json').write_text(json.dumps(load_provider(r),indent=2,allow_nan=False));return p
