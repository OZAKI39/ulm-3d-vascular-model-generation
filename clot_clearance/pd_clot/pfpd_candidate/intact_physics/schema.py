"""Explicit method/material/unit identity; rejected inputs never trigger a solve."""
import copy,hashlib,json,math,re
from pathlib import Path
import numpy as np
from ..curved_clot_audit.common import ROOT,OLD,read,digest
METHODS=('FEM_NATIVE_P2','FEM_NATIVE_P3_CHECK','PD2_VARIATIONAL_SVK_SAVED')
MATERIAL='SVK_DEMO_1000PA_LAMBDA0'
USE='EXPLORATORY_QUASISTATIC_INTACT_DEMO'
class RequestError(ValueError):
    def __init__(self,code,detail):super().__init__(detail);self.code=code
def reject(code,detail):raise RequestError(code,detail)
def finite_tree(x):
    if isinstance(x,dict):
        for k,v in x.items():
            key=k.lower()
            if any(word in key for word in ['frequency','pulse_width','wall_velocity','time_series','time_sequence']) or key in ['t','times','t_sequence','dt','frequency_hz','velocity','freq_hz']:
                reject('DYNAMIC_MODEL_NOT_IMPLEMENTED',f'{k} is not a quasistatic load coefficient')
            finite_tree(v)
    elif isinstance(x,(list,tuple)):
        for v in x:finite_tree(v)
    elif isinstance(x,(int,float)) and not isinstance(x,bool):
        if not math.isfinite(x):reject('INVALID_NUMERIC_INPUT','Nonfinite numeric request')
def validate_request(request):
    if isinstance(request,dict) and request.get("schema_version")=="3.0":
        from .v3 import validate_request as validate_v3
        return validate_v3(request)
    if isinstance(request,dict) and request.get("schema_version")=="2.0":
        from .v2 import validate_request as validate_v2
        return validate_v2(request)
    r=copy.deepcopy(request)
    if not isinstance(r,dict):reject('INVALID_REQUEST','JSON object required')
    finite_tree(r)
    required=['schema_version','request_id','method_id','geometry_id','geometry','material_id','material','units','wall_mode','wall','load_mode','load','response_spec','solver','output_dir','shared_inputs','task_configuration','intended_use','allowed_error','action','source']
    if not set(required)<=set(r) or set(r)-set(required)-{'accuracy_tier'}:reject('INVALID_REQUEST',f'Required explicit keys differ: missing={set(required)-set(r)}, unknown={set(r)-set(required)}')
    if r.get('accuracy_tier','TIER_A') not in ['TIER_A','TIER_B','TIER_C']:reject('INVALID_USE_REQUIREMENT','Unknown accuracy tier')
    if r['schema_version']!='1.0':reject('UNSUPPORTED_SCHEMA','Only 1.0')
    if not isinstance(r['request_id'],str) or not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]{0,79}',r['request_id']):reject('INVALID_REQUEST_ID','Use a simple unique case ID')
    if r['method_id'] not in METHODS:reject('UNSUPPORTED_FOR_THIS_INTERFACE_VERSION','Unknown method; no fallback backend')
    if r['units']!=dict(length='m',pressure='Pa',displacement='m',energy='J'):reject('INVALID_UNITS','This version requires SI m, Pa, m, J')
    if r['material_id']!=MATERIAL or r['material']!=dict(mu_Pa=1000.,lambda_Pa=0.,provenance='DEMONSTRATION_PARAMETERS_NOT_BIOCALIBRATED'):reject('UNSUPPORTED_MATERIAL','Only the exact uncalibrated demo SVK card is executable')
    C=read(r['task_configuration'])
    if r['geometry_id']!=C['geometry_id'] or r['geometry']!=dict(mesh_file=C['mesh_file'],sha256=C['mesh_sha256'],R_m=C['R_m'],h1_m=C['h1_m']):reject('GEOMETRY_IDENTITY_MISMATCH','Only the explicitly identified accepted fine P1 geometry is supported')
    if Path(C['mesh_file']).resolve()!=(OLD/'geometry/fine/mesh.npz').resolve() or digest(C['mesh_file'])!=C['mesh_sha256']:reject('GEOMETRY_IDENTITY_MISMATCH','Accepted mesh bytes differ')
    if r['wall_mode'] not in ['FIXED_DISK','PRESCRIBED_INPLANE_STRETCH']:reject('UNSUPPORTED_WALL','No adhesion, contact or dynamic wall law')
    if set(r['wall'])!={'epsilon_w','origin_m'} or r['wall']['origin_m']!=[0.,0.,0.]:reject('INVALID_WALL','Reference coordinate origin must remain fixed')
    e=r['wall']['epsilon_w']
    if isinstance(e,bool) or not isinstance(e,(float,int)) or abs(e)>.01:reject('INVALID_WALL','This small-deformation exploratory interface bounds |epsilon_w| <= 0.01; not biological validity')
    if r['wall_mode']=='FIXED_DISK' and e!=0:reject('INVALID_WALL','Fixed disk requires epsilon_w=0')
    if r['load_mode'] not in ['GLOBAL','PATCH']:reject('UNSUPPORTED_LOAD','Only GLOBAL/PATCH nominal dead traction')
    if set(r['load'])!={'p_Pa','tau_x_Pa','tau_y_Pa','center_m','ell_m'}:reject('INVALID_LOAD','Unexpected load parameter or missing mode coefficient')
    for key in ['p_Pa','tau_x_Pa','tau_y_Pa','ell_m']:
        if isinstance(r['load'][key],bool) or not isinstance(r['load'][key],(int,float)):reject('INVALID_LOAD','Scalar SI mode coefficients/width required')
    center=np.asarray(r['load']['center_m'],float);R=C['R_m']
    if center.shape!=(3,) or not np.isfinite(center).all() or center[2]<0 or abs(np.linalg.norm(center)-R)>1e-10*R:reject('INVALID_PATCH_CENTER','Center must lie on ideal exposed hemisphere in reference coordinates')
    if not 0<r['load']['ell_m']<=2*R:reject('INVALID_PATCH_WIDTH','Require 0 < ell <= 2R in this interface')
    if not isinstance(r['intended_use'],str) or not r['intended_use'].strip():reject('INVALID_USE_REQUIREMENT','Explicit nonempty intended_use required')
    if r['intended_use']!=USE and r['action']!='REUSE_RESULT':reject('UNSUPPORTED_USE','New solve use remains restricted; read-only reuse may record unapproved new use context')
    if set(r['allowed_error'])!={'absolute','relative','source'}:reject('INVALID_USE_REQUIREMENT','Use explicit optional tolerance/source fields')
    for key in ['absolute','relative']:
        x=r['allowed_error'][key]
        if x is not None and (not isinstance(x,(float,int)) or isinstance(x,bool) or x<0):reject('INVALID_USE_REQUIREMENT','A supplied tolerance must be a nonnegative number')
    if any(r['allowed_error'][k] is not None for k in ['absolute','relative']) and not r['allowed_error']['source']:reject('INVALID_USE_REQUIREMENT','A tolerance needs a real stated source; this module does not approve it')
    expected=dict(sampling_id=C['sampling']['sampling_id'],measure_id=C['sampling']['measure_id'],profiles=['axis','transverse'],native_quadrature_extrema=True)
    if r['response_spec']!=expected:reject('UNSUPPORTED_RESPONSE_SPEC','Use the frozen measures/profiles; no new samples or diagnostics by default')
    if r['solver']!=C['solver']:reject('UNSUPPORTED_SOLVER_CONFIGURATION','Use frozen solver configuration; no silent adaptation')
    if r['action'] not in ['SOLVE','LOAD_SAVED_BASELINE','REUSE_RESULT']:reject('UNSUPPORTED_ACTION','Unsupported action')
    if r['method_id']=='PD2_VARIATIONAL_SVK_SAVED' and (r['action']=='SOLVE' or r['load_mode']!='GLOBAL' or r['wall_mode']!='FIXED_DISK' or [r['load'][k] for k in ['p_Pa','tau_x_Pa','tau_y_Pa']]!=[3.,3.,0.]):reject('UNSUPPORTED_FOR_THIS_INTERFACE_VERSION','PD2 local-load/new-material/new-wall solving is not implemented; never relabel FEM as PD2')
    if r['action']=='SOLVE':
        expected_cases={'LP':(2,'PATCH',3.,0.,'FIXED_DISK',0.),'LS':(2,'PATCH',0.,3.,'FIXED_DISK',0.),'LM':(2,'PATCH',3.,3.,'FIXED_DISK',0.),'W':(2,'GLOBAL',0.,0.,'PRESCRIBED_INPLANE_STRETCH',1e-3),'LM_P3':(3,'PATCH',3.,3.,'FIXED_DISK',0.)}
        if r['request_id'] not in expected_cases:reject('OUTSIDE_FROZEN_CAMPAIGN_BUDGET','Only five declared scientific cases may be solved this round')
        degree=3 if r['method_id']=='FEM_NATIVE_P3_CHECK' else 2
        values=(degree,r['load_mode'],r['load']['p_Pa'],r['load']['tau_x_Pa'],r['wall_mode'],e)
        if values!=expected_cases[r['request_id']] or r['load']['tau_y_Pa']!=0 or r['load']['center_m']!=C['patch']['center_m'] or r['load']['ell_m']!=C['patch']['ell_m']:reject('OUTSIDE_FROZEN_CAMPAIGN_BUDGET','No changed geometry/material/load/wall campaign is authorized')
    if r['action']!='SOLVE' and not r['source']:reject('MISSING_SOURCE','Saved result/state source required')
    out=Path(r['output_dir']).resolve();shared=Path(r['shared_inputs']).resolve()
    if out.name!=r['request_id'] or out.parent!=shared.parent/'cases':reject('INVALID_OUTPUT_PATH','Output must be unique task cases/request_id path')
    if shared!=Path(r['task_configuration']).resolve().parent/'shared':reject('INVALID_SHARED_INPUT_PATH','Shared inputs must belong to this task')
    return r
def physical_identity(r):return {k:r[k] for k in ['method_id','geometry_id','geometry','material_id','material','wall_mode','wall','load_mode','load','response_spec']}
def identity_hash(r):return hashlib.sha256(json.dumps(physical_identity(r),sort_keys=True).encode()).hexdigest()
