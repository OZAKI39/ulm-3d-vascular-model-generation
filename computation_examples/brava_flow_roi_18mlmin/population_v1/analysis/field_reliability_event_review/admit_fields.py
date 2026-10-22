from common import *
import datetime
import subprocess
import pyvista as pv

def triangle_set(poly):
    xyz=np.asarray(poly.points,dtype=np.float64)
    points,inv=np.unique(xyz,axis=0,return_inverse=True)
    faces=np.sort(inv[poly.faces.reshape(-1,4)[:,1:]],axis=1)
    faces=faces[np.lexsort(faces.T[::-1])]
    return points,faces

def main():
    assert (HERE/'EVENT_FIELD_PROTOCOL.json').exists()
    inputs=read(HERE/'PROTOCOL_SOURCE_SHA256.json')
    def reg(p):
        inputs[str(p.relative_to(BASE))]=sha(p);return inputs[str(p.relative_to(BASE))]
    names=['WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03']
    original={k:triangle_set(pv.read(paths('M0')[0]/f'SV_MESH/mesh-surfaces/{k}.vtp')) for k in names}
    original_exec=read(paths('M0')[1]/'reports/execution.json')
    fields={}
    for name in NAMES:
        mesh,case=paths(name)
        meta={}
        for rel in ['run/solver.xml','run/PETSC_OPTIONS.txt','reports/execution.json','reports/progress.json','policy.json']:
            if (case/rel).exists():reg(case/rel)
        for p in [meshpath(name),flowpath(name)]: reg(p)
        assert flowpath(name).stat().st_size>1000000
        roles={}
        for k in names:
            path=mesh/f'SV_MESH/mesh-surfaces/{k}.vtp';reg(path)
            current=triangle_set(pv.read(path)); old=original[k]
            identical=all(np.array_equal(a,b) for a,b in zip(old,current))
            roles[k]=dict(exact_coordinate_and_triangle_sets=identical,points=len(current[0]),triangles=len(current[1]),sha256=inputs[str(path.relative_to(BASE))])
            assert identical,(name,k,'GEOMETRY_CONFOUNDED')
        execution=read(case/'reports/execution.json'); progress=read(case/'reports/progress.json')
        policy=read(case/'policy.json')
        actual_mesh=inputs[str(meshpath(name).relative_to(BASE))];actual_flow=inputs[str(flowpath(name).relative_to(BASE))]
        if name=='M0':
            assert actual_mesh=='e755a6ce04eca518fc02df5cd6139c98887e94617bdf9410d72f55b99f443c94'
            assert actual_flow=='73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc'
        else:
            ident=read(case/'STUDY_IDENTITY.json');reg(case/'STUDY_IDENTITY.json')
            assert ident['mesh_sha256']==actual_mesh and ident['flow_sha256']==actual_flow
            qa=mesh.parent/'reports/FIXED_PLC_AUDIT.json';reg(qa);q=read(qa)
            gate=mesh.parent/'reports/GATE_MESH.json';reg(gate)
            meta['geometry_QA_status']=q['status'];meta['raw_validity']=q['raw']
            assert q['status'].startswith('PASS') and not q['failures']
        solver_equal=(case/'run/solver.xml').read_bytes()==(paths('M0')[1]/'run/solver.xml').read_bytes()
        assert solver_equal
        assert execution['exit_code']==0 and execution['status']=='PASS'
        assert execution['solver_sha256']==original_exec['solver_sha256']
        assert execution['PETSc_library_sha256']==original_exec['PETSc_library_sha256']
        assert execution['initial_state']=='zero'
        assert len(progress['intervals'])>=5
        assert all(r['E_u']<=policy['velocity_change_limit'] and r['E_Q']<=policy['flow_change_limit'] for r in progress['intervals'][-5:])
        native=case/execution['final_vtu'];reg(native)
        checkpoint=Path(execution['checkpoint']['path'])
        if checkpoint.exists():reg(checkpoint)
        history=execution.get('history',{})
        hist_summary={k:v for k,v in history.items() if not isinstance(v,(list,dict))}
        fields[name]=dict(status='ADMITTED_METADATA_AND_COMPLETE_BOUNDARY',
          geometry=dict(roles=roles,full_boundary_exact=True,coordinate_transform='identity; retained coordinates in SI original frame',**meta),
          BC=dict(exact_original_solver_XML=solver_equal,profile='Flat; Impose_flux=true; Zero_out_perimeter=true; Q outward=-3e-7 m3/s',
            outlets_Pa=[556.97065554884477,0,17.137427153562019],outlet_semantics='Neu natural traction in same executable',wall='Dir zero'),
          physics=dict(equation='steady limit of incompressible Newtonian Navier-Stokes',rho=1056,mu=0.00345312),
          solver_discretization=dict(solver_sha256=execution['solver_sha256'],PETSc_library_sha256=execution['PETSc_library_sha256'],
             same_P1_stabilization_executable_and_XML=True,dt_s=.01,spectral_radius=.5,MPI_ranks=execution['MPI_ranks'],
             PETSC_OPTIONS=execution['PETSC_OPTIONS'],same_options=execution['PETSC_OPTIONS']==original_exec['PETSC_OPTIONS']),
          steady=dict(status=execution['status'],exit_code=execution['exit_code'],final_step=execution['final_step'],
             stopping_record=progress['stop'],last_five_intervals=progress['intervals'][-5:],
             original_thresholds={k:policy[k] for k in ['steady_last_intervals','velocity_change_limit','flow_change_limit','mass_limit']},
             solver_history_summary=hist_summary,local_weak_residual='NOT_EVALUABLE: saved solver aggregate residual only; no spatial stabilized weak residual reconstructed'),
          raw_payload=dict(mesh=str(meshpath(name).relative_to(BASE)),flow=str(flowpath(name).relative_to(BASE)),mesh_sha256=actual_mesh,flow_sha256=actual_flow,
             native=str(native.relative_to(BASE)),native_sha256=inputs[str(native.relative_to(BASE))],initial_state='zero',actual_independent_solve=True,
             native_mapping='pending sampling-stage actual array check'),valid_event_probes='pending',
          remaining_confounders=['Same P1 method is not independent truth; mesh family is locally refined and not uniform asymptotic sequence.','Domain geometry fixed; changes in interior mesh over full domain are not isolated replacement of one cell.'])
        print('ADMITTED',name,actual_mesh,actual_flow,flush=True)
    # Register required original authoritative records and bounded prior evidence.
    for rel in ['reports/ACTIVE_FLOW.json','population_v1/runs/realization_0000_contactfix/IDENTITY.json',
      'population_v1/analysis/speed_event_cause_review/CAUSE_RESULTS.json','population_v1/analysis/speed_event_cause_review/INPUT_SHA256.json',
      'population_v1/analysis/speed_event_cause_review/EXECUTION_RECEIPT.json','population_v1/analysis/speed_event_cause_review/REPRODUCE.md',
      'population_v1/analysis/acceleration_diagnosis/WALL_RELEASE_AND_FEM.json','population_v1/analysis/failure_speed_review_20260929/largest_increment_context.json',
      'population_v1/cfd_wall_audit/WALL_CFD_AUDIT_REPORT_ZH.md','mesh_particle_validation_v2/FINAL_P1_SUMMARY.json']:
        reg(BASE/rel)
    put('INPUT_SHA256.json',inputs)
    put('FIELD_COMPARABILITY.json',dict(fields=fields,admission_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
      new_CFD_count=0,reused_CFD_count=3,reused_CFD=NAMES,reference_payload_count=2,
      optional_CFD='NOT_TRIGGERED: both existing fixed-PLC fields available, geometry/BC/steady admitted.'))
    resources={k:(Path('/sys/fs/cgroup')/k).read_text().strip() for k in ['cpu.max','cpuset.cpus.effective','memory.max','memory.current','memory.stat','memory.events']}
    resources['GPU']=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.total,memory.used,memory.free,utilization.gpu','--format=csv'],text=True)
    resources['disk']=subprocess.check_output(['df','-h',str(BASE)],text=True)
    resources['time_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    put('RESOURCES.json',resources)

if __name__=='__main__':main()
