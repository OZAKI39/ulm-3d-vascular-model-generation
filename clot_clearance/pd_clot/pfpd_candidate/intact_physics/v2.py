"""Explicit new-version entry point; does not bypass or relax v1 solve guards."""
import time,traceback
from pathlib import Path
from .schema import RequestError,finite_tree
from .api import export_result,enforce_quality_contract
from ..curved_clot_audit.common import read,dump,digest
from ..materials_literature.ogden import CARD

def validate_request(r):
    finite_tree(r)
    required={'schema_version','request_id','method_id','material_id','material','geometry','load','wall','intended_use','allowed_error','task_configuration','output_dir','baseline_task','action'}
    if r.get('action') in ['RESUME_LINEAR_SOLVER_RECOVERY','RESUME_IMPLEMENTATION_BUG']:required.add('recovery')
    if set(r)!=required or r['schema_version']!='2.0':raise RequestError('INVALID_V2_REQUEST','Explicit v2 schema keys required')
    if r['material_id']!=CARD['material_id'] or r['material']!=CARD:raise RequestError('MATERIAL_IDENTITY_MISMATCH','Only the source-identified Ogden card is authorized')
    C=read(r['task_configuration']);O=Path(r['task_configuration']).resolve().parent
    if O!=Path(__file__).resolve().parents[3]/'verification/clot_literature_material_and_controlled_fracture_v1':raise RequestError('CAMPAIGN_IDENTITY_MISMATCH','Wrong task path')
    expected_out=O/'cases'/r['request_id']
    if r['action']=='RESUME_LINEAR_SOLVER_RECOVERY':
        rec=read(O/'LINEAR_SOLVER_RECOVERY.json')
        if r['request_id']!='MAT_W' or r['recovery']!=rec or digest(rec['checkpoint'])!=rec['checkpoint_sha256']:raise RequestError('INVALID_RECOVERY','Only recorded same-case trusted checkpoint')
        expected_out=Path(rec['output_dir'])
    if r['action']=='RESUME_IMPLEMENTATION_BUG':
        rec=read(O/'FRACTURE_IMPLEMENTATION_RECOVERY.json')
        if r['request_id']!='FRACTURE_PATH' or r['recovery']!=rec or digest(rec['last_accepted_checkpoint']['path'])!=rec['last_accepted_checkpoint']['sha256']:raise RequestError('INVALID_RECOVERY','Only recorded last trusted fracture checkpoint')
    if Path(r['output_dir']).resolve()!=expected_out:raise RequestError('INVALID_OUTPUT_PATH','Use new task cases and case ID')
    if r['action'] not in ['SOLVE','RESUME_LINEAR_SOLVER_RECOVERY','RESUME_IMPLEMENTATION_BUG']:raise RequestError('UNSUPPORTED_V2_ACTION','This v2 execution adapter only solves authorized new cases; old reuse remains v1')
    if r['request_id'] not in ['MAT_LM','MAT_W','FRACTURE_PATH']:raise RequestError('OUTSIDE_V2_BUDGET','Only the two hemisphere states and one fracture path')
    if r['request_id'].startswith('MAT_'):
        if r['method_id']!=C['hemisphere']['method_id']:raise RequestError('METHOD_IDENTITY_MISMATCH','No implicit backend change')
        old=read(Path(r['baseline_task'])/'requests/LM.json')
        if r['geometry']!=old['geometry'] or digest(r['geometry']['mesh_file'])!=r['geometry']['sha256']:raise RequestError('GEOMETRY_IDENTITY_MISMATCH','Reuse exact accepted fine mesh')
        expected=(3.,3.,0.,0.) if r['request_id']=='MAT_LM' else (0.,0.,0.,.1)
        expected_load=old['load'].copy();expected_load.update(p_Pa=expected[0],tau_x_Pa=expected[1],tau_y_Pa=expected[2])
        if r['load']!=expected_load or r['wall']!=dict(epsilon_w=expected[3]):raise RequestError('OUTSIDE_V2_BUDGET','No extra or silently ignored load/wall fields')
        if (r['load']['p_Pa'],r['load']['tau_x_Pa'],r['load']['tau_y_Pa'],r['wall']['epsilon_w'])!=expected or r['load']['center_m']!=old['load']['center_m'] or r['load']['ell_m']!=old['load']['ell_m']:raise RequestError('OUTSIDE_V2_BUDGET','Changed new material loading not authorized')
    else:
        G=read(O/'FRACTURE_GEOMETRY_CONFIGURATION.json')
        if r['load']!=C['fracture'] or r['wall']!=dict(grip=G['grip'],symmetry=G['symmetry']) or r['geometry']!=dict(mesh_file=str(O/'shared/specimen/final_mesh.npz'),sha256=digest(O/'shared/specimen/final_mesh.npz')):raise RequestError('FRACTURE_CONFIGURATION_MISMATCH','Only frozen path/mesh/BC permitted')
        if r['method_id']!=C['fracture']['method_id'] or not read(O/'FRACTURE_READINESS.json')['ready']:raise RequestError('FRACTURE_NOT_RUN_BLOCKED','Fracture requires actual passing readiness')
    return r

def run(r):
    result=dict(schema_version='2.0',request_id=r.get('request_id'),method_id=r.get('method_id'),material_id=r.get('material_id'),success=False,response=None,solve_status='NOT_STARTED',intended_use=r.get('intended_use'),use_requirement=r.get('allowed_error'),use_approval='NOT_BIOLOGICALLY_VALIDATED',downstream_control_tolerances='NEEDS_USE_REQUIREMENT');out=None;owned=False
    try:
        validate_request(r);out=Path(r['output_dir'])
        if r['action']=='RESUME_IMPLEMENTATION_BUG':
            if (out/'result.json').exists() or (out/'RESUME_STARTED_001').exists():raise RequestError('DUPLICATE_REQUEST_ID','Recovery cannot overwrite a completed or started recovery')
            owned=True;(out/'RESUME_STARTED_001').write_text(str(time.time()));dump(out/'resume_request_001.json',r)
        else:
            if out.exists():raise RequestError('DUPLICATE_REQUEST_ID','Successful or started case cannot be overwritten')
            out.mkdir(parents=True);owned=True;(out/'STARTED').write_text(str(time.time()));dump(out/'request.json',r)
        if r['request_id'].startswith('MAT_'):
            from ..materials_literature.hemisphere import run_case
        else:
            from ..fem_fracture_candidate.path import run_case
        result.update(run_case(r));result=enforce_quality_contract(result);return export_result(result,out)
    except Exception as exc:
        result.update(success=False,response=None,solve_status=getattr(exc,'code','IMPLEMENTATION_OR_NUMERICAL_FAILURE'),error=dict(type=type(exc).__name__,message=str(exc)))
        if not owned:out=Path(r['task_configuration']).parent/'rejected_requests'/str(time.time_ns());out.mkdir(parents=True,exist_ok=False)
        (out/'failure.log').write_text(traceback.format_exc());return export_result(result,out)
