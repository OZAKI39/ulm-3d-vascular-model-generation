"""One gated V2 manufactured solve adapter and raw production P1 field audit."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import sys

import numpy as np
import pyvista as pv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NAME = 'R1_GMSH_V2'
CASE = ROOT/'manufactured_pipe'/NAME/'cfd'
sys.path[:0] = [str(ROOT/'scripts'),str(ROOT.parent/'vendor')]
from pipe_measurements import PipeMeasurements
from pipe_probe_geometry import frozen_probes, strict_cells
from particle_field_inputs import curl, shear, errors
from poiseuille import reference


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def put(path, value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def prepare():
    from prepare_pipe_cfd import main as prepare_existing
    qa=read(ROOT/'GMSH_REFERENCE_MESH_V2_RESULT.json')
    assert qa['status']=='GMSH_REFERENCE_MESH_V2_PASS' and all(qa['checks'].values())
    plan=read(HERE/'GMSH_REFERENCE_MESH_V2_PLAN.json')
    mesh=HERE/'candidate_v2/mesh'
    assert sha(mesh/'volume.vtu')==qa['volume_sha256']
    proto=read(ROOT/'manufactured_pipe/PROTOCOL.json')
    exact=reference(np.array([[0.,0.,0.]]),plan['radius_m'],plan['length_m'],proto['mean_velocity_m_s'],proto['mu_Pa_s'])
    inputs=dict(name='manufactured/'+NAME,radius_m=plan['radius_m'],length_m=plan['length_m'],
                h_bulk_m=plan['bulk_size_m'],policy_name='GMSH_V2_VALIDATION_ONLY',
                analytic=dict(Q_m3_s=exact['Q_m3_s'],WSS_Pa=exact['WSS_Pa'],pressure_gradient_Pa_m=exact['pressure_gradient_Pa_m']),
                model_sha256=qa['identity']['cad_sha256'],mesh_sha256=qa['volume_sha256'],
                plan_sha256=sha(HERE/'GMSH_REFERENCE_MESH_V2_PLAN.json'),
                volume_mesher='Frozen Gmsh 4.15.2 HXT, first-order tetra',
                surface_design=plan['surface_meshing'])
    inp=ROOT/'straight_pipe/manufactured'/NAME
    inp.mkdir(exist_ok=False)
    put(inp/'INPUT.json',inputs)
    # Existing adapter needs this receipt. Reuse the actual V2 finite-facet/cap
    # tests, without a second geometry study or any mesh modification.
    put(mesh/'FIDELITY_QA.json',dict(status='PASS',surface=dict(PASS=qa['checks']['finite_facet_geometry'],
        finite_facet_max_deviation_m=qa['central_finite_facet_deviation_max_m']),
        caps=qa['cap_analytical_trace_diagnostic'],source='../../../GMSH_REFERENCE_MESH_V2_RESULT.json',
        source_sha256=sha(ROOT/'GMSH_REFERENCE_MESH_V2_RESULT.json'),method='Actual V2 triangles; exact finite-facet bound and degree-4 cap quadrature'))
    prepare_existing(NAME,mesh,manufactured=True)
    policy=read(CASE/'policy.json')
    policy['maximum_wall_time_s']=1800
    put(CASE/'policy.json',policy)
    resource=read(CASE/'reports/RESOURCE_PLAN.json')
    resource.update(hard_wall_budget_s=1800,hard_tree_RSS_GiB=24,
        runtime_authorization='User explicitly authorizes one solve with 1800 s external hard limit. Original conservative estimate is retained; it is not a runtime promise.')
    put(CASE/'reports/RESOURCE_PLAN.json',resource)
    hashes=read(CASE/'input_hashes.json')
    hashes['policy.json']=sha(CASE/'policy.json')
    put(CASE/'input_hashes.json',hashes)
    # Numerical settings must equal the existing manufactured solver case.
    import xml.etree.ElementTree as ET
    new=ET.parse(CASE/'run/solver.xml')
    old=ET.parse(ROOT/'manufactured_pipe/R1_P0/cfd/run/solver.xml')
    assert ET.tostring(new.getroot())==ET.tostring(old.getroot())
    assert (CASE/'run/PETSC_OPTIONS.txt').read_bytes()==(ROOT/'manufactured_pipe/R1_P0/cfd/run/PETSC_OPTIONS.txt').read_bytes()
    put(CASE/'reports/FROZEN_SOLVER_IDENTITY.json',dict(solver_xml_sha256=sha(CASE/'run/solver.xml'),
        PETSC_OPTIONS_sha256=sha(CASE/'run/PETSC_OPTIONS.txt'),manufactured_R1_P0_numerical_settings_equal=True,
        rho_kg_m3=proto['rho_kg_m3'],mu_Pa_s=proto['mu_Pa_s'],Umean_m_s=proto['mean_velocity_m_s'],
        MPI_ranks=8,OMP_NUM_THREADS=1,hard_limit_s=1800))
    print('PREPARED_ONE_MANUFACTURED_CFD',flush=True)


def analyze(case=CASE, supervision_name='manufactured_gmsh_v2'):
    CASE = case
    qa=read(ROOT/'GMSH_REFERENCE_MESH_V2_RESULT.json')
    assert qa['status']=='GMSH_REFERENCE_MESH_V2_PASS'
    supervision=read(ROOT/'monitoring'/supervision_name/'EXECUTION.json')
    if not (CASE/'reports/execution.json').exists():
        assert supervision['stop_reason'] is not None or supervision['exit_code'] != 0
        diagnostic=read(CASE/'reports/SERIALIZED_OUTPUT_MONITOR.json')
        history=(CASE/'run/8-procs/histor.dat').read_text()
        timing=[dict(step=int(m[1]),iteration=int(m[2]),cumulative_native_s=float(m[3]))
                for m in re.finditer(r'^\s*NS\s+(\d+)-(\d+)\w*\s+([\d.eE+-]+)',history,re.MULTILINE)]
        previous=0.
        for row in timing:
            row['increment_native_s']=row['cumulative_native_s']-previous;previous=row['cumulative_native_s']
        flux=(CASE/'run/8-procs/B_NS_Velocity_flux.txt').read_text()
        completed_steps=[int(m[1]) for m in re.finditer(r'^\s*(\d+)\s+[\d.eE+-]+\s',flux,re.MULTILINE)]
        put(ROOT/'FIELD_VALIDATION.json',dict(status='REFERENCE_PIPE_FIELD_FAIL',
            reason='The single CFD did not produce a completed production convergence receipt',
            failure_stage='CFD_RUNTIME_OR_EXECUTION',stop_reason=supervision['stop_reason'],
            supervised_exit_code=supervision['exit_code'],elapsed_s=supervision['elapsed_s'],
            supervision_sha256=sha(ROOT/'monitoring/manufactured_gmsh_v2/EXECUTION.json'),
            peak_tree_RSS_GiB=supervision['peak_tree_RSS_GiB'],native_timing=timing,
            last_completed_native_step=max(completed_steps,default=0),
            timing_note='Cumulative T and derived increments are limited to the native histor.dat printed precision; wall clock guard is authoritative.',
            serialization_monitor=diagnostic,required_consecutive_steady_intervals=5,
            source_sha256={str(p.relative_to(ROOT)):sha(p) for p in [CASE/'run/solver.log',CASE/'run/8-procs/histor.dat',
                CASE/'reports/SERIALIZED_OUTPUT_MONITOR.json',CASE/'reports/FROZEN_SOLVER_IDENTITY.json',ROOT/'scripts/pipe_measurements.py',Path(__file__)]},
            mesh_status=qa['status'],field_metrics=None,CFD_launches=1,
            local_field_accuracy_assessed=False,field_certified=False,particle_cases_executed=0,
            automatic_BraVa=False))
        return
    execution=read(CASE/'reports/execution.json')
    if execution['status']!='PASS':
        put(ROOT/'FIELD_VALIDATION.json',dict(status='REFERENCE_PIPE_FIELD_FAIL',reason='Production convergence gate failed',
            failure_stage='CFD_CONVERGENCE',solver_status=execution['status'],
            solver_execution_sha256=sha(CASE/'reports/execution.json'),
            local_field_accuracy_assessed=False,field_certified=False))
        return
    proto=read(ROOT/'manufactured_pipe/PROTOCOL.json')
    inp=read(ROOT/'straight_pipe/manufactured'/NAME/'INPUT.json')
    radius,length=inp['radius_m'],inp['length_m']
    mean,mu=proto['mean_velocity_m_s'],proto['mu_Pa_s']
    mesh=CASE/'SV_MESH/mesh-complete.mesh.vtu'
    assert sha(mesh)==qa['volume_sha256']
    field=CASE/execution['final_vtu']
    measure=PipeMeasurements(CASE/'SV_MESH/mesh_arrays.npz',inp['analytic']['Q_m3_s'],mean)
    u,p=measure.read(field)
    points,tetra=measure.points,measure.tetra
    grid=pv.read(mesh)
    dx=points[tetra[:,1:]]-points[tetra[:,:1]]
    G=np.linalg.solve(dx,u[tetra[:,1:]]-u[tetra[:,:1]]).transpose(0,2,1)
    levels,x,groups,planes,angles=frozen_probes(proto,radius,length)
    with np.load(HERE/'candidate_v2/QA_SAMPLES.npz') as frozen:
        assert np.array_equal(x,frozen['probes_m'])
    cells,locator=strict_cells(grid,points,tetra,x)
    valid=cells>=0
    c=cells[valid]
    b=np.linalg.solve(dx[c].transpose(0,2,1),(x[valid]-points[tetra[c,0]])[...,None])[...,0]
    weights=np.c_[1-b.sum(1),b]
    uh=np.full((len(x),3),np.nan);gh=np.full((len(x),3,3),np.nan)
    uh[valid]=np.einsum('ni,nij->nj',weights,u[tetra[c]]);gh[valid]=G[c]
    exact=reference(x,radius,length,mean,mu)
    n=-np.c_[np.cos(angles),np.sin(angles),np.zeros(len(x))]
    actual=dict(u=uh,G=gh,s=shear(gh,n),curl=curl(gh))
    expected=dict(u=exact['u_m_s'],G=exact['G_s_inv'],s=shear(exact['G_s_inv'],n),curl=curl(exact['G_s_inv']))
    # Diagnostic decomposition only: no substitution for the numerical field.
    nodal=reference(points,radius,length,mean,mu)['u_m_s']
    nodal[measure.wall_nodes]=0.
    iu=np.full_like(uh,np.nan)
    iu[valid]=np.einsum('ni,nij->nj',weights,nodal[tetra[c]])
    gi=np.linalg.solve(dx,nodal[tetra[:,1:]]-nodal[tetra[:,:1]]).transpose(0,2,1)
    rows=[]
    for group,d in enumerate(levels):
        mask=groups==group;role=proto['probe_roles'][group];complete=bool(valid[mask].all())
        metrics={k:errors(actual[k][mask],expected[k][mask]) if complete else None for k in actual}
        gates=({k:complete and metrics[k]['relative_L2']<=limit for k,limit in [('u',.05),('G',.08),('s',.08),('curl',.08)]}
               if role=='PARTICLE_RELEVANT' else dict(u=complete and metrics['u']['relative_L2']<=.05) if role=='BULK' else {})
        rows.append(dict(wall_distance_m=float(d),role=role,samples=int(mask.sum()),inside=int(valid[mask].sum()),
            metrics=metrics,gates=gates,
            analytical_P1_interpolant_u_L2=errors(iu[mask],expected['u'][mask])['relative_L2'] if complete else None,
            CFD_minus_interpolant_u_L2=errors(uh[mask]-iu[mask]+expected['u'][mask],expected['u'][mask])['relative_L2'] if complete else None,
            analytical_P1_interpolant_G_L2=errors(gi[cells[mask]],expected['G'][mask])['relative_L2'] if complete else None))
    grid.point_data['audit_u']=u;grid.point_data['audit_p']=p
    sections=[]
    for z in proto['analysis_z_over_L']:
        section=grid.slice(normal=[0,0,1],origin=[0,0,z*length]).triangulate()
        tri=section.faces.reshape(-1,4)[:,1:];xyz=section.points[tri]
        area=np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)/2
        q=float(np.dot(area,section['audit_u'][tri,2].mean(1)))
        pp=float(np.dot(area,section['audit_p'][tri].mean(1))/area.sum())
        sections.append(dict(z_over_L=z,z_m=z*length,Q_m3_s=q,Q_relative_error=abs(q/inp['analytic']['Q_m3_s']-1),pressure_mean_Pa=pp))
    slope=float(np.polyfit([s['z_m'] for s in sections],[s['pressure_mean_Pa'] for s in sections],1)[0])
    dp_error=abs(slope/inp['analytic']['pressure_gradient_Pa_m']-1)
    with np.load(HERE/'candidate_v2/mesh/BOUNDARY_TOPOLOGY.npz') as topology:
        mask=topology['facet_tags']==1
        wall=points[topology['boundary_triangles'][mask]];owners=topology['adjacent_tetra'][mask]
    normals=np.cross(wall[:,1]-wall[:,0],wall[:,2]-wall[:,0]);area=np.linalg.norm(normals,axis=1)/2
    normals/=2*area[:,None]
    traction=np.einsum('nij,nj->ni',mu*(G[owners]+G[owners].transpose(0,2,1)),normals)
    wss=np.linalg.norm(traction-normals*np.sum(traction*normals,axis=1)[:,None],axis=1)
    centers=wall.mean(1);central=(centers[:,2]>=.4*length)&(centers[:,2]<=.6*length)
    wss_error=float(np.sqrt(np.dot(area[central],(wss[central]-exact['WSS_Pa'])**2)/area[central].sum())/exact['WSS_Pa'])
    gates=dict(solver=True,complete_frozen_probes=bool(valid.all()),finite=bool(np.isfinite(u).all() and np.isfinite(p).all() and np.isfinite(G).all()),
        particle_bands=all(all(r['gates'].values()) for r in rows if r['role']=='PARTICLE_RELEVANT'),
        bulk_velocity=all(r['gates']['u'] for r in rows if r['role']=='BULK'),
        internal_Q=all(s['Q_relative_error']<=.05 for s in sections),pressure_gradient=dp_error<=.08)
    out=CASE/'analysis';out.mkdir(exist_ok=False)
    np.savez_compressed(out/'FIELD_SAMPLES.npz',points_m=x,cell=cells,groups=groups,planes=planes,
        velocity=uh,gradient=gh,exact_velocity=expected['u'],exact_gradient=expected['G'],
        interpolant_velocity=iu,gradient_all_tetra=G,nodal_velocity=u,nodal_pressure=p)
    sources=[mesh,field,CASE/'reports/execution.json',CASE/'reports/FROZEN_SOLVER_IDENTITY.json',
             HERE/'GMSH_REFERENCE_MESH_V2_PLAN.json',ROOT/'PARTICLE_CFD_INPUT_CONTRACT.json',out/'FIELD_SAMPLES.npz',Path(__file__)]
    result=dict(status='REFERENCE_PIPE_FIELD_PASS' if all(gates.values()) else 'REFERENCE_PIPE_FIELD_FAIL',VALIDATION_ONLY=True,
        case=NAME,checks=gates,groups=rows,internal_sections=sections,fitted_dp_dz_Pa_m=slope,dp_dz_relative_error=dp_error,
        central_area_weighted_raw_WSS_relative_L2=wss_error,WSS_role='diagnostic only',
        source_sha256={str(f.relative_to(ROOT)):sha(f) for f in sources},
        raw_P1_G_definition='G_ij=du_i/dx_j; no smoothing/recovery; inward exact-cylinder n for both fields; s=(I-nn^T)Gn',
        error_definition='Per-band Euclidean vector/Frobenius relative L2 over all 96 frozen probes; no zero-component percentages',
        interpolation_diagnostic='u_h-u_exact=(u_h-Iu_exact)+(Iu_exact-u_exact); norms do not add. Interpolant never replaces CFD.',
        field_frozen=all(gates.values()),particle_oracle_launched=False,automatic_BraVa=False)
    put(ROOT/'FIELD_VALIDATION.json',result)
    print(json.dumps(dict(status=result['status'],checks=gates,particle_bands=[r for r in rows if r['role']=='PARTICLE_RELEVANT']),indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','analyze'])
    {'prepare':prepare,'analyze':analyze}[parser.parse_args().stage]()
