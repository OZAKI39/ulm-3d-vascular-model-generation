"""Independent physical verification and lossless new frozen-flow export."""
from pathlib import Path
import sys, json, hashlib, shutil, csv, xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.sv13 import stop_gate
from sv_validation.sv13n import checkpoint_one_rank
from sv_validation.sv11 import linear_gate, nonlinear_gate
CASE=ROOT/'flow_cases/mean-2p0-mmps'
Q=1.551359160885543e-14

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,data): path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')

def check_bc(path,baseline):
    tree=ET.parse(path); old=ET.parse(baseline)
    inlet=tree.find(".//Add_BC[@name='INLET']")
    assert float(inlet.findtext('Value'))==-Q, 'Signed INLET flux is wrong'
    assert {k:inlet.findtext(k) for k in ('Type','Profile','Impose_flux','Zero_out_perimeter')}==dict(
        Type='Dir',Profile='Flat',Impose_flux='true',Zero_out_perimeter='true')
    for name in ('OUTLET_01','OUTLET_02','OUTLET_03','WALL'):
        bc=tree.find(f".//Add_BC[@name='{name}']")
        assert bc.findtext('Type')==('Dir' if name=='WALL' else 'Neu')
        assert float(bc.findtext('Value'))==0
    inlet.find('Value').text=old.find(".//Add_BC[@name='INLET']/Value").text
    assert ET.tostring(tree.getroot())==ET.tostring(old.getroot()), 'Unexpected solver/physics/mesh configuration change'
    return True

def check_mass(m):
    incoming=-m['signed_outward_boundary_flows_m3_s']['INLET']
    out=m['outlet_flows_m3_s']
    assert set(out)=={'OUTLET_01','OUTLET_02','OUTLET_03'}
    assert np.isfinite([incoming,*out.values()]).all()
    assert incoming>0 and all(v>0 for v in out.values()), 'Inlet or outlet direction is wrong'
    assert abs(incoming-Q)/Q<=1e-6, 'Target inflow mismatch'
    assert abs(sum(out.values())-incoming)/Q<=1e-6, 'Flux conservation failure'
    return True

def mesh_matches(grid,measure):
    assert np.array_equal(grid.points,measure.points), 'Point coordinates/order changed'
    cells=grid.cells.reshape(-1,5)
    assert np.all(cells[:,0]==4) and np.all(grid.celltypes==10), 'Non-tetrahedral cells'
    assert np.array_equal(np.sort(cells[:,1:],axis=1),np.sort(measure.tetra,axis=1)), 'Tetrahedra changed'
    return True

def run():
    import pyvista as pv
    reports=CASE/'reports';baseline=ROOT/'frozen_reference'
    policy=json.loads((CASE/'policy.json').read_text())
    execution=json.loads((reports/'execution.json').read_text())
    check_bc(CASE/'run/solver.xml',baseline/'run/solver.xml')
    for name,expected in json.loads((reports/'old_baseline_hashes.json').read_text()).items():
        assert sha(baseline/name)==expected, 'Old frozen baseline changed: '+name
    for name,expected in json.loads((CASE/'input_hashes.json').read_text()).items():
        assert sha(CASE/name)==expected, 'Case input changed: '+name
    assert execution['status']=='PASS'
    assert sha(CASE/'run/solver.log')==execution['log_sha256'], 'Native log changed'
    linear_gate(execution); nonlinear_gate(execution['history'])
    final=CASE/execution['final_vtu'];measure=SolutionMeasurements(CASE/'SV_MESH/mesh_arrays.npz',Q,.002)
    u,p=measure.read(final); grid=pv.read(final);mesh_matches(grid,measure)
    m=measure.measure(u,p);check_mass(m)
    assert m['wall_noslip_pass'] and m['velocity_finite'] and m['pressure_finite']
    assert stop_gate(execution['intervals'],m,0,0,policy)
    qualified=execution['stop']['step']
    prior=final.with_name('result_%03d.vtu'%qualified)
    # Native VTU suffix can have a different padding width; find numerically.
    prior=next(x for x in final.parent.glob('result_*.vtu') if int(x.stem.rsplit('_',1)[1])==qualified)
    uq,pq=measure.read(prior);mq=measure.measure(uq,pq)
    final_change=measure.velocity_l2(u-uq)/measure.velocity_l2(u)
    final_q_change=max(abs(m['signed_outward_boundary_flows_m3_s'][r]-mq['signed_outward_boundary_flows_m3_s'][r])/Q for r in m['signed_outward_boundary_flows_m3_s'])
    assert final_change<=1e-5 and final_q_change<=1e-6
    cp=final.with_name('stFile_%03d.bin'%execution['final_step'])
    checkpoint=checkpoint_one_rank(cp,execution['final_step'],policy['dt_s'])
    speed=np.linalg.norm(u,axis=1)
    # Four-point, degree-two tetra quadrature applied to |P1 vector field|.
    # Magnitude is not polynomial: clearly report this representative mean as quadrature.
    weights=np.full((4,4),.1381966011250105);np.fill_diagonal(weights,.5854101966249685)
    local_u=np.einsum('qi,tij->tqj',weights,u[measure.tetra])
    mean_speed=float(np.dot(measure.volumes,np.linalg.norm(local_u,axis=2).mean(axis=1))/measure.volumes.sum())
    m.update(velocity_min_m_s=float(speed.min()),velocity_volume_mean_m_s=mean_speed,
        velocity_nodal_mean_m_s=float(speed.mean()),velocity_volume_rms_m_s=measure.velocity_l2(u)/np.sqrt(measure.volumes.sum()),
        representative_mean_method='4-point tetra quadrature of norm(P1 velocity); SI volume weights',
        source='actual_vtu_surface_integration',inlet_area_m2=policy['A_in_m2'],
        inlet_actual_mean_m_s=m['Q_in_m3_s']/policy['A_in_m2'])
    old_u,old_p=measure.read(baseline/'flow/steady_flow_stage_sv1_3q.vtu')
    old=measure.measure(old_u,old_p)
    ratio=m['Q_in_m3_s']/old['Q_in_m3_s']
    comparison=dict(old_Q_in_m3_s=old['Q_in_m3_s'],flow_ratio=ratio,
        old_max_speed_mm_s=old['velocity_max_m_s']*1000,new_max_speed_mm_s=m['velocity_max_m_s']*1000,
        maximum_speed_ratio=m['velocity_max_m_s']/old['velocity_max_m_s'],
        velocity_relative_L2_difference_from_scaled_old=measure.velocity_l2(u-ratio*old_u)/measure.velocity_l2(u),
        old_outlet_fractions=old['outlet_fractions'],new_outlet_fractions=m['outlet_fractions'],
        note='Scaled old field is used only for a comparison metric, never as the solved or exported field.')
    frozen=CASE/'frozen_flow'; frozen.mkdir(exist_ok=True)
    shutil.copyfile(final,frozen/'steady_flow_mean_2p0_mmps.vtu')
    shutil.copyfile(cp,frozen/cp.name)
    np.savez_compressed(frozen/'flow_arrays_si.npz',points_m=measure.points,tetra=measure.tetra,
        boundary_triangles=measure.boundary,facet_tags=measure.tags,velocity_m_s=u,pressure_pa=p)
    exported=pv.read(frozen/'steady_flow_mean_2p0_mmps.vtu');mesh_matches(exported,measure)
    assert sha(final)==sha(frozen/'steady_flow_mean_2p0_mmps.vtu')
    assert np.array_equal(exported['Velocity'],u) and np.array_equal(exported['Pressure'].reshape(-1),p)
    with np.load(frozen/'flow_arrays_si.npz') as a:
        assert np.array_equal(a['velocity_m_s'],u) and np.array_equal(a['pressure_pa'],p)
        assert np.array_equal(a['tetra'],measure.tetra)
    manifest=dict(case='mean-2p0-mmps',status='PASS',raw_solver_vtu=str(final.relative_to(CASE)),
        Q_target_m3_s=Q,inlet_mean_velocity_m_s=.002,density_kg_m3=1056.,dynamic_viscosity_pa_s=.00345312,
        mesh_sha256=sha(CASE/'SV_MESH/mesh-complete.mesh.vtu'),node_count=len(measure.points),tetra_count=len(measure.tetra),
        final_step=execution['final_step'],physical_time_s=execution['final_step']*policy['dt_s'],
        units=dict(points='m',velocity='m/s',pressure='Pa',flow='m3/s'),
        mesh_note='NPZ uses the original positive-oriented tetra order. VTU preserves solver-local vertex order; same point order and same tetra sets were checked.',
        arrays=dict(points='points_m',tetra='tetra',velocity='velocity_m_s',pressure='pressure_pa',boundaries='boundary_triangles',tags='facet_tags'),
        face_ids={'WALL':1,'OUTLET_03':2,'OUTLET_01':3,'INLET':4,'OUTLET_02':5},
        files={x.name:dict(sha256=sha(x),bytes=x.stat().st_size) for x in sorted(frozen.iterdir()) if x.is_file() and x.name!='manifest.json'},
        source_solver_sha256=execution['solver_sha256'],candidate_for_particle_background=True,
        rigid_wall_newtonian=True,experimental_validation=False)
    dump(frozen/'manifest.json',manifest)
    result=dict(status='PASS',boundary_sanity=True,old_baseline_unchanged=True,case_inputs_unchanged=True,
        measurements=m,comparison=comparison,steady=dict(first_qualifying_step=qualified,
        final_step=execution['final_step'],intervals=execution['intervals'],
        final_relative_velocity_change=final_change,final_normalized_flow_change=final_q_change),
        native_checkpoint=checkpoint,policy=policy,export=manifest)
    dump(reports/'physics_validation.json',result)
    with (reports/'boundary_flows.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['boundary','signed_outward_m3_s','signed_outward_pL_s','outlet_fraction'])
        for role,value in m['signed_outward_boundary_flows_m3_s'].items():
            writer.writerow([role,format(value,'.17g'),format(value*1e15,'.17g'),m['outlet_fractions'].get(role,'')])
    with (reports/'steady_intervals.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=execution['intervals'][0].keys());writer.writeheader();writer.writerows(execution['intervals'])
    print(json.dumps(dict(status='PASS',measurements=m,comparison=comparison),indent=2))

if __name__=='__main__':run()
