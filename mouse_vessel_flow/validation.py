"""Read-only evidence checks and small network regressions; no solver/trajectory launch."""
import ast,importlib,json,sys,traceback
from pathlib import Path
import numpy as np
from .paths import PROJECT_ROOT,CASE_NAMES,project_path,case_config,new_output
from .provenance import sha256,protected_originals,verify_source_integrity

def independent_imports():
    modules=[];failures=[]
    for path in sorted((PROJECT_ROOT/'mouse_vessel_flow').rglob('*.py')):
        name='.'.join(path.relative_to(PROJECT_ROOT).with_suffix('').parts)
        if name.endswith('.__init__'):name=name[:-9]
        try:importlib.import_module(name);modules.append(name)
        except Exception as exc:failures.append(dict(module=name,error=repr(exc)))
    outside=[]
    for name,module in tuple(sys.modules.items()):
        if name.startswith('mouse_vessel_flow'):
            f=getattr(module,'__file__',None)
            if f and not Path(f).resolve().is_relative_to(PROJECT_ROOT):outside.append(f)
    links=[str(p.relative_to(PROJECT_ROOT)) for p in PROJECT_ROOT.rglob('*')
           if '.venv' not in p.parts and '.git' not in p.parts and p.is_symlink()]
    return dict(status='FAIL' if failures or outside or links else 'PASS',
                modules=len(modules),failures=failures,external_project_modules=outside,project_symlinks=links,
                working_directory=str(Path.cwd()),python=sys.executable)

def copied_files():
    # Original migration evidence remains immutable; this manifest tracks the
    # retained dataset after explicitly authorized removals.
    rows=json.loads(project_path('configs/active_data_manifest.json').read_text());fail=[]
    for row in rows:
        p=project_path(row['path'])
        if not p.is_file() or sha256(p)!=row['sha256']:fail.append(row['path'])
    return dict(status='FAIL' if fail else 'PASS',files=len(rows),bytes=sum(r['size'] for r in rows),failures=fail)

def field_check(name):
    from .flow_field import load_case
    from .fem import inspect_case
    from .particles.audit import validate_pair
    env=load_case(name,with_particles=False)
    tetra,_=validate_pair(env.mesh,env.flow)
    # P1 interpolation at element centroids is the mean of the four original nodal values.
    cells=np.linspace(0,len(tetra)-1,48,dtype=int);points=np.asarray(env.mesh.points,dtype=np.float64)[tetra[cells]].mean(axis=1)
    maximum=0.;gradient_error=0.
    for cell,x in zip(cells,points):
        sample=env.field.sample(x);expected=np.asarray(env.flow.point_data['Velocity'])[tetra[cell]].mean(axis=0)
        np.testing.assert_allclose(sample.velocity_m_s,expected,rtol=2e-10,atol=1e-15)
        assert sample.inside_lumen
        maximum=max(maximum,float(np.max(abs(sample.velocity_m_s-expected))))
        gradient_error=max(gradient_error,float(np.max(abs(sample.velocity_gradient_s_inv-env.field.gradients_s_inv[cell]))))
    cfg=case_config(name);base=project_path(cfg['data_root'])/'flow'
    acceptance={}
    if (base/'reports/execution.json').exists():
        saved=json.loads((base/'reports/execution.json').read_text());expected=saved['final']
        assert saved['status']=='PASS'
        for port,q in expected['outlet_flows_m3_s'].items():
            np.testing.assert_allclose(env.audit[port]['signed_Q_m3_s'],q,rtol=2e-14,atol=0)
        np.testing.assert_allclose(env.audit['INLET']['signed_Q_m3_s'],expected['Q_in_m3_s'],rtol=2e-14,atol=0)
        np.testing.assert_array_equal([env.field.pressure_nodes_pa.min(),env.field.pressure_nodes_pa.max()],expected['pressure_range_pa'])
        cp=base/'run/1-procs/stFile_last.bin';cp_info=saved['checkpoint']
        acceptance=dict(status='PASS',historical_physical_acceptance=saved['status'],outlet_fluxes_match=True,
            pressure_range_matches=True,stop=saved['stop'])
    elif name in ('mean_2mmps','stage_q'):
        cpbase=project_path(cfg['fem_reference'])/'frozen_reference' if name=='mean_2mmps' else base
        cp_info=json.loads((cpbase/'fem_checkpoint/checkpoint_manifest.json').read_text())
        cp=cpbase/'fem_checkpoint/stFile_071.bin'
    else:raise ValueError('Unrecognized checkpoint provenance')
    from .solver_support.checkpoint_checks import checkpoint_one_rank
    assert sha256(cp)==cp_info['sha256']
    checkpoint=checkpoint_one_rank(cp,cp_info['step'],inspect_case(name)['dt_s'])
    return dict(status='PASS',sha256=env.new_flow_sha256,points=env.mesh.n_points,cells=env.mesh.n_cells,
        flux=env.audit,centroid_samples=len(cells),maximum_interpolation_abs_error_m_s=maximum,
        maximum_gradient_abs_error_s_inv=gradient_error,solver=inspect_case(name),
        native_checkpoint=checkpoint,acceptance=acceptance,
        pressure_range_pa=[float(env.field.pressure_nodes_pa.min()),float(env.field.pressure_nodes_pa.max())],
        velocity_max_m_s=float(np.linalg.norm(env.field.velocity_nodes_m_s,axis=1).max()))

def wss_check():
    import pyvista as pv
    from .solver_support.wss_case import recover
    cfg=case_config('roi_only');flow=pv.read(project_path(cfg['flow']['path']))
    with np.load(project_path('inputs/reference/fem/frozen_reference/SV_MESH/mesh_arrays.npz')) as z:
        current,_=recover(z['points_m'],z['tetra'],z['boundary_triangles'],z['facet_tags'],
                          np.asarray(flow['Velocity']),np.asarray(flow['Pressure']),cfg['dynamic_viscosity_pa_s'])
    saved=pv.read(project_path('results/roi_only/visualization/input_data/field_diagnostics/data/wall_wss_si.vtp'))
    errors={}
    for key in ('WSS_raw_Pa','Tangential_viscous_traction_Pa','Outward_normal','Area_m2','Parent_tetra_zero_based'):
        a=np.asarray(current.cell_data[key]);b=np.asarray(saved.cell_data[key])
        np.testing.assert_allclose(a,b,rtol=5e-13,atol=1e-25)
        errors[key]=float(np.max(abs(a-b)))
    return dict(status='PASS',wall_facets=saved.n_cells,max_abs_differences=errors,
        coordinates='Original float64 mesh_arrays.npz; no geometry replacement',
        original_recorded_validation='validation/wss_validation_v2/WSS_VALIDATION_REPORT_V2.md')

def validate(case,output,*,full_integrity=False):
    out=new_output(output);checks={}
    def check(key,fn):
        try:checks[key]=fn()
        except Exception as exc:checks[key]=dict(status='FAIL',error=repr(exc),traceback=traceback.format_exc())
        print(key+': '+checks[key]['status'],file=sys.stderr,flush=True)
    check('imports_and_paths',independent_imports)
    check('runtime_source_identity',lambda:dict(status='PASS',count=len(verify_source_integrity())))
    check('original_protected_source_identity',lambda:dict(status='PASS',count=len(protected_originals())))
    if full_integrity:check('all_copied_data',copied_files)
    from .network import validate_design
    check('roi_network_and_h0_regression',validate_design)
    check('wss_frozen_recovery',wss_check)
    from .particles.analysis import verify_tracks
    for name in (CASE_NAMES if case=='all' else (case,)):
        check(name+'_field',lambda name=name:field_check(name))
        if case_config(name).get('particles'):check(name+'_particles',lambda name=name:verify_tracks(name))
    checks['new_formal_CFD']=dict(status='NOT_RUN',reason='Frozen workflow restructuring; no new solve authorized')
    from .native_environment import native_status
    check('native_CFD_environment',lambda:native_status('native_cfd'))
    check('native_mesher_environment',lambda:native_status('native_mesher'))
    checks['new_formal_particle_batch']=dict(status='NOT_RUN',reason='Original frozen trajectories are verified read-only')
    checks['fine_vessel_and_O2_sensitivity']=dict(status='NOT_RUN',original_status='SKIPPED_BY_USER',reason='Cancellation preserved')
    report=dict(status='FAIL' if any(c['status']=='FAIL' for c in checks.values()) else 'PASS',
        scope='LOCAL_REFACTOR_VALIDATION_ONLY',checks=checks,
        limits='PASS applies only to executed checks; native compute availability and new simulations are reported separately.')
    (out/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
    return report
