"""Versioned, six-case Level-1 screening entry point; older guards unchanged."""
import copy
import hashlib
import json
import re
import time
import traceback
from pathlib import Path
from .schema import RequestError, finite_tree
from .policy import decorate_result, validate_metadata, POLICY_ID
from ..curved_clot_audit.common import ROOT, read, dump, digest
from ..materials_literature.ogden import CARD

TASK = ROOT/'verification/policy_and_fast_material_screening_v1'
CASES = {
    'G_P': ('GLOBAL',3.,0.,0.), 'G_S': ('GLOBAL',0.,3.,0.),
    'L_P': ('PATCH',3.,0.,0.), 'L_S': ('PATCH',0.,3.,0.),
    'L_M': ('PATCH',3.,3.,0.), 'WALL': ('GLOBAL',0.,0.,.01),
}

def physical_identity(r):
    return {k:r[k] for k in ['method_id','geometry_id','geometry','material_id','material',
                             'wall_mode','wall','load_mode','load','response_spec']}

def identity_hash(r):
    return hashlib.sha256(json.dumps(physical_identity(r),sort_keys=True).encode()).hexdigest()

def validate_request(request):
    if not isinstance(request,dict):
        raise RequestError('INVALID_V3_REQUEST','JSON object required')
    r=copy.deepcopy(request);finite_tree(r)
    required={'schema_version','request_id','method_id','geometry_id','geometry','material_id',
              'material','units','wall_mode','wall','load_mode','load','response_spec','solver',
              'output_dir','shared_inputs','task_configuration','baseline_task','intended_use',
              'allowed_error','accuracy_tier','action','source'}
    if set(r)!=required or r['schema_version']!='3.0':
        raise RequestError('INVALID_V3_REQUEST','Explicit schema 3.0 fields required')
    if not isinstance(r['request_id'],str) or not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]{0,79}',r['request_id']):
        raise RequestError('INVALID_REQUEST_ID','Simple task-local request ID required')
    if r['action'] not in ['SOLVE','REUSE_RESULT']:
        raise RequestError('UNSUPPORTED_V3_ACTION','No dynamic/fracture/PD2 action')
    if r['method_id']!='FEM_NATIVE_P2':
        raise RequestError('METHOD_IDENTITY_MISMATCH','Only native P2 FEM screening')
    if r['material_id']!=CARD['material_id'] or r['material']!=CARD:
        raise RequestError('MATERIAL_IDENTITY_MISMATCH','Unchanged identified Ogden parameters required')
    if Path(r['task_configuration']).resolve()!=TASK/'CONFIGURATION.json':
        raise RequestError('CAMPAIGN_IDENTITY_MISMATCH','New task configuration required')
    C=read(r['task_configuration'])
    if digest(r['task_configuration'])!=read(TASK/'FREEZE_IDENTITY.json')['configuration_sha256']:
        raise RequestError('CONFIGURATION_CHANGED','Frozen six-case configuration changed')
    if r['geometry_id']!=C['geometry_id'] or r['geometry']!=C['geometry'] or digest(r['geometry']['mesh_file'])!=r['geometry']['sha256']:
        raise RequestError('GEOMETRY_IDENTITY_MISMATCH','Accepted fine mesh bytes required')
    if r['units']!=dict(length='m',pressure='Pa',displacement='m',energy='J'):
        raise RequestError('INVALID_UNITS','SI units required')
    if r['solver']!=C['solver'] or r['response_spec']!=C['response_spec']:
        raise RequestError('UNSUPPORTED_CONFIGURATION','Frozen solver and measures required')
    if Path(r['shared_inputs']).resolve()!=TASK/'shared' or Path(r['baseline_task']).resolve()!=Path(C['baseline_task']):
        raise RequestError('INPUT_IDENTITY_MISMATCH','Frozen common measure and historical source required')
    if not isinstance(r['intended_use'],str) or not r['intended_use'].strip() or r['accuracy_tier'] not in ['TIER_A','TIER_B','TIER_C']:
        raise RequestError('INVALID_USE_REQUIREMENT','Explicit use and tier required')
    ae=r['allowed_error']
    if not isinstance(ae,dict) or set(ae)!={'absolute','relative','source'}:
        raise RequestError('INVALID_USE_REQUIREMENT','Explicit tolerance/source required')
    for key in ['absolute','relative']:
        v=ae[key]
        if v is not None and (isinstance(v,bool) or not isinstance(v,(float,int)) or v<0):
            raise RequestError('INVALID_USE_REQUIREMENT','Nonnegative tolerance required')
    if not ae['source']:
        raise RequestError('INVALID_USE_REQUIREMENT','Policy source required')
    if Path(r['output_dir']).resolve()!=TASK/'cases'/r['request_id']:
        raise RequestError('INVALID_OUTPUT_PATH','Task-owned case output required')
    if r['action']=='SOLVE':
        if r['request_id'] not in CASES:
            raise RequestError('OUTSIDE_V3_BUDGET','Exactly six declared application cases')
        mode,p,t,e=CASES[r['request_id']]
        if r['load_mode']!=mode or r['load']!=dict(p_Pa=p,tau_x_Pa=t,tau_y_Pa=0.,center_m=C['patch']['center_m'],ell_m=C['patch']['ell_m']) or r['wall']!=dict(epsilon_w=e,origin_m=[0.,0.,0.]) or r['wall_mode']!=('PRESCRIBED_INPLANE_STRETCH' if e else 'FIXED_DISK'):
            raise RequestError('OUTSIDE_V3_BUDGET','Changed load or wall parameters are outside this campaign')
        if r['accuracy_tier']!='TIER_A' or ae!=C['allowed_error']:
            raise RequestError('UNAPPROVED_SOLVE_USE','This campaign authorizes Tier-A screening only')
        decision=read(TASK/'MATERIAL_CONVENTION_DECISION.json')
        if not decision['authorized'] or not read(TASK/'MATERIAL_READINESS.json')['ready']:
            raise RequestError('MATERIAL_NOT_READY','Source, material and native variational gates required')
    elif not r['source']:
        raise RequestError('MISSING_SOURCE','Successful source result required')
    return r

def run(request):
    from .api import enforce_quality_contract, export_result
    start=time.perf_counter();r=request if isinstance(request,dict) else {};owned=False;out=None
    result=dict(schema_version='3.0',request_id=r.get('request_id'),method_id=r.get('method_id'),
                material_id=r.get('material_id'),geometry_id=r.get('geometry_id'),load_id=r.get('request_id'),
                success=False,response=None,solve_status='NOT_STARTED',new_equilibrium_solve_calls=0,
                quantity_definition='SI; reference-volume energy and RMS; current-position force/moment',
                reference_scope='SINGLE_FINE_GEOMETRY_P2; NO_NEW_BIOLOGICAL_OR_PD2_VALIDATION')
    try:
        t=time.perf_counter();r=validate_request(r);validation_s=time.perf_counter()-t
        out=Path(r['output_dir']);out.mkdir(parents=True,exist_ok=False);owned=True
        (out/'STARTED').open('x').write(str(time.time()));dump(out/'request.json',r)
        result.update(physical_identity=physical_identity(r),identity_sha256=identity_hash(r))
        if r['action']=='REUSE_RESULT':
            original=read(r['source'])
            if not original.get('success') or original.get('identity_sha256')!=identity_hash(r):
                raise RequestError('RESULT_REUSE_IDENTITY_MISMATCH','Successful matching physical state required')
            for a in original.get('artifacts',{}).values():
                if digest(a['path'])!=a['sha256']:
                    raise RequestError('RESULT_REUSE_IDENTITY_MISMATCH','Source artifact changed')
            result=copy.deepcopy(original)
            result.update(request_id=r['request_id'],solve_status='REUSED_WITHOUT_SOLVE',
                          source_use_requirement=dict(intended_use=original.get('intended_use'),allowed_error=original.get('use_requirement'),source=original.get('source'),accuracy_tier=original.get('accuracy_tier')),
                          use_approval='NOT_GRANTED_BY_REUSE',new_equilibrium_solve_calls=0,
                          reused_result_source=dict(path=r['source'],sha256=digest(r['source'])))
        else:
            from .screening_backend import prepare, solve_or_load, evaluate
            p=prepare(r);p['profile']['request_input_validation_s']=validation_s
            solution=solve_or_load(p);result.update(evaluate(p,solution))
        enforce_quality_contract(result);decorate_result(result,r);validate_metadata(result)
        result['interface_wall_seconds']=time.perf_counter()-start
        return export_result(result,out)
    except Exception as exc:
        result.update(success=False,response=None,solve_status=getattr(exc,'code','IMPLEMENTATION_OR_NUMERICAL_FAILURE'),error=dict(type=type(exc).__name__,message=str(exc)))
        if not owned:
            out=TASK/'rejected_requests'/str(time.time_ns());out.mkdir(parents=True,exist_ok=False)
        result['new_equilibrium_solve_calls']=int((out/'EQUILIBRIUM_CALLED').exists())
        decorate_result(result,r);result['interface_wall_seconds']=time.perf_counter()-start
        (out/'failure.log').write_text(traceback.format_exc())
        return export_result(result,out)
