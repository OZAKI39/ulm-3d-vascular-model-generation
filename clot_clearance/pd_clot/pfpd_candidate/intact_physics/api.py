"""Callable and serializable interface; backend identity never changes implicitly."""
import copy,json,time,traceback,resource,re
from pathlib import Path
import numpy as np
from .schema import validate_request,RequestError,identity_hash,physical_identity,USE
from ..curved_clot_audit.common import read,dump,digest,OLD
def prepare_inputs(request):
    if request.get('schema_version')=='3.0':
        from .screening_backend import prepare
        return prepare(request)
    from .fem_backend import prepare
    return prepare(request)
def solve_or_load(prepared):
    if isinstance(prepared,dict) and prepared.get('schema_version')=='3.0':
        from .screening_backend import solve_or_load as execute
        return execute(prepared)
    from .fem_backend import solve_or_load as execute
    return execute(prepared)
def evaluate_response(prepared,solution):
    if isinstance(prepared,dict) and prepared.get('schema_version')=='3.0':
        from .screening_backend import evaluate
        return evaluate(prepared,solution)
    from .fem_backend import evaluate
    return evaluate(prepared,solution)
def export_result(result,output_dir):
    dump(Path(output_dir)/'result.json',result);return read(Path(output_dir)/'result.json')
def metadata(r):
    r={k:v if isinstance(v,str) else None for k,v in r.items()}
    return dict(schema_version='1.0',request_id=r.get('request_id'),method_id=r.get('method_id'),geometry_id=r.get('geometry_id'),material_id=r.get('material_id'),wall_id=r.get('wall_mode'),load_id=r.get('load_mode'),quantity_definition=None,reference_scope='NOT_ASSESSED',discrepancy='historical PD2/FEM output-specific differences retained; no new cross-method agreement claim',solve_status='NOT_STARTED',data_source_identity=None,use_label=USE,downstream_control_tolerances='NEEDS_USE_REQUIREMENT',Ftilde=dict(status='NOT_APPLICABLE',value=None),success=False,response=None)
def enforce_quality_contract(result):
    """Keep failed candidate numbers diagnostic, never an approved response."""
    if not result.get('success') or result.get('solve_status')=='NUMERICAL_QUALITY_FAILED' or result.get('quality') is False:
        result['success']=False
        if result.get('response') is not None:result['diagnostic_response']=result.pop('response')
        result['response']=None
        result.setdefault('error',dict(type='NumericalQualityFailure',message='One or more required numerical gates did not pass',failed_gates=result.get('failed_gates',['NUMERICAL_QUALITY_FAILED_UNSPECIFIED_GATE'])))
    return result
def run(request):
    if isinstance(request,dict) and request.get("schema_version")=="3.0":
        from .v3 import run as run_v3
        return run_v3(request)
    if isinstance(request,dict) and request.get("schema_version")=="2.0":
        from .v2 import run as run_v2
        return run_v2(request)
    start=time.perf_counter(); result=metadata(request if isinstance(request,dict) else {});out=None;owned=False
    try:
        r=validate_request(request);out=Path(r['output_dir'])
        if out.exists():raise RequestError('DUPLICATE_REQUEST_ID','Refusing to overwrite an existing case directory')
        out.mkdir(parents=True,exist_ok=False);owned=True;dump(out/'request.json',r);(out/'STARTED').open('x').write(time.ctime()+'\n')
        result.update(physical_identity=physical_identity(r),identity_sha256=identity_hash(r),input_coefficients=r['load'],wall_parameters=r['wall'],request_action=r['action'],use_requirement=r['allowed_error'])
        if r['action']=='REUSE_RESULT':
            original=read(r['source'])
            if original.get('identity_sha256')!=identity_hash(r) or not original.get('success'):raise RequestError('RESULT_REUSE_IDENTITY_MISMATCH','Reuse requires matching physical identity and an actual successful source result')
            for a in original.get('artifacts',{}).values():
                if not Path(a['path']).is_file() or digest(a['path'])!=a['sha256']:raise RequestError('RESULT_REUSE_IDENTITY_MISMATCH','Source artifact changed or missing')
            result=copy.deepcopy(original)
            result.update(source_use_requirement=dict(intended_use=original.get('intended_use',original.get('use_label')),allowed_error=original.get('use_requirement'),source=original.get('source'),accuracy_tier=original.get('accuracy_tier')),intended_use=r['intended_use'],use_label=r['intended_use'],use_requirement=copy.deepcopy(r['allowed_error']),source=r['source'],use_approval='NOT_GRANTED_BY_REUSE',downstream_control_tolerances='NEEDS_USE_REQUIREMENT_ASSESSMENT',request_id=r['request_id'],request_action='REUSE_RESULT',solve_status='REUSED_WITHOUT_SOLVE',reused_result_source=dict(path=r['source'],sha256=digest(r['source'])),new_equilibrium_solve_calls=0)
        elif r['method_id']=='PD2_VARIATIONAL_SVK_SAVED':
            expected=OLD/'cases/fine_M/state.npz'
            if Path(r['source']).resolve()!=expected.resolve():raise RequestError('STATE_IDENTITY_MISMATCH','This saved PD2 adapter supports the original fine M state only')
            z=np.load(expected);s=read(expected.parent/'summary.json');assert np.array_equal(z['X'],np.load(r['geometry']['mesh_file'])['X'])
            result.update(success=bool(s['PASS']),solve_status='SAVED_BASELINE_LOADED',data_source_identity=dict(path=str(expected),sha256=digest(expected)),quantity_definition='historical PD2 own nodal-volume RMS, original apex and energy; no new common-volume field evaluation',reference_scope='SAVED_PD2_FINE_M_ONLY; historical comparisons unchanged',response=dict(apex_displacement_m=s['apex_displacement_m'],nodal_u_RMS_m=s['u_RMS_m'],energy_J=s['energy_J'],balance=s['balance']),Ftilde=dict(status='SAVED_BOND_ASSOCIATED_DEFINITION_NOT_REEVALUATED',value=None),artifacts=dict(saved_state=dict(path=str(expected),sha256=digest(expected))),new_equilibrium_solve_calls=0)
        else:
            prepared=prepare_inputs(r);solution=solve_or_load(prepared);response=evaluate_response(prepared,solution);result.update(response)
        result=enforce_quality_contract(result)
        from .policy import decorate_result
        decorate_result(result,r)
        result['interface_wall_seconds']=time.perf_counter()-start;result['interface_peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        (out/'interface.log').write_text(json.dumps(dict(request_id=r['request_id'],method_id=r['method_id'],status=result['solve_status'],wall_seconds=result['interface_wall_seconds']))+'\n')
        return export_result(result,out)
    except Exception as exc:
        result.update(success=False,solve_status=getattr(exc,'code','IMPLEMENTATION_OR_NUMERICAL_FAILURE'),error=dict(type=type(exc).__name__,message=str(exc)),response=None,interface_wall_seconds=time.perf_counter()-start)
        # Never replace a previous result; preflight rejections get their own small record.
        if not owned or (out/'result.json').exists():
            cfg=request.get('task_configuration') if isinstance(request,dict) else None
            parent=Path(cfg).resolve().parent if isinstance(cfg,(str,Path)) and Path(cfg).is_file() else Path.cwd()
            safe_id=re.sub(r'[^A-Za-z0-9_.-]','_',str(result.get('request_id') or 'invalid'))[:80]
            out=parent/'rejected_requests'/(str(time.time_ns())+'_'+safe_id);out.mkdir(parents=True,exist_ok=False)
        else:out.mkdir(parents=True,exist_ok=True)
        result['new_equilibrium_solve_calls']=int((out/'EQUILIBRIUM_CALLED').exists())
        from .policy import decorate_result
        decorate_result(result,request)
        (out/'failure.log').write_text(traceback.format_exc());dump(out/'result.json',result)
        return result
