"""Derive strain/shear rates from the identified final FEM solution in SI units."""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json
import numpy as np
import pyvista as pv

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def load_core():
    path=HERE/'audit/wss_core.py'
    expected=json.loads((ROOT/'input_data/field_diagnostics/COMPUTE_VALIDATION.json').read_text())['core_sha256']
    assert sha(path)==expected,'Use the exact production P1 gradient implementation'
    spec=importlib.util.spec_from_file_location('production_wss',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def rates(gradient):
    d=.5*(gradient+np.swapaxes(gradient,-1,-2))
    norm=np.sqrt(np.sum(d*d,axis=(-2,-1)))
    shear=np.sqrt(2*np.sum(d*d,axis=(-2,-1)))
    return d,norm,shear

def analytical_checks(core):
    # Exact affine velocities on a physical tetrahedron exercise the same P1
    # gradient and scalar operations as the actual field, including SI scaling.
    x=np.array([[0,0,0],[2,0,0],[0,3,0],[0,0,4]],float)*1e-6
    t=np.array([[0,1,2,3]])
    cases=[('simple_shear',np.array([[0,1000,0],[0,0,0],[0,0,0]],float),1000.),
           ('rigid_rotation',np.array([[0,-700,0],[700,0,0],[0,0,0]],float),0.),
           ('incompressible_extension',np.diag([800.,-400.,-400.]),np.sqrt(3)*800.)]
    results=[]
    for name,g,expected in cases:
        u=x@g.T+np.array([.001,.002,-.003])
        recovered=core.p1_gradients(x,t,u)
        d,norm,actual=rates(recovered)
        error=abs(float(actual[0])-expected)
        assert error<1e-9 and np.max(abs(recovered[0]-g))<1e-9
        results.append(dict(case=name,expected_s_inv=expected,actual_s_inv=float(actual[0]),
                            absolute_error_s_inv=error,passed=True))
    return results

def summary(values,weight,region,quantity):
    order=np.argsort(values);cdf=(np.cumsum(weight[order])-.5*weight[order])/weight.sum()
    qs=np.interp([.05,.5,.95],cdf,values[order])
    return dict(region=region,quantity=quantity,unit='s^-1',count=len(values),
                weight='tetra_volume' if region=='whole_volume' else 'triangle_area',
                min=float(values.min()),mean=float(np.average(values,weights=weight)),
                P5=float(qs[0]),P50=float(qs[1]),P95=float(qs[2]),max=float(values.max()))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE)
    args=parser.parse_args();out=args.output.resolve();(out/'data').mkdir(parents=True,exist_ok=True)
    if (out/'COMPUTE_VALIDATION.json').exists():raise FileExistsError('Preserve completed derived data')
    lock=json.loads((HERE/'audit/source_lock.json').read_text())
    for name,digest in lock.items():assert sha(ROOT/name)==digest,name
    compute=json.loads((ROOT/'input_data/field_diagnostics/COMPUTE_VALIDATION.json').read_text())
    core=load_core();checks=analytical_checks(core)
    source=ROOT/'input_data/frozen_flow/flow_arrays_si.npz'
    field=ROOT/'input_data/frozen_flow/steady_flow_mean_2p0_mmps.vtu'
    assert sha(source)==compute['flow_arrays_sha256']
    assert sha(field)==compute['source_field_sha256']
    a=np.load(source);x=a['points_m'];tetra=a['tetra'];tri=a['boundary_triangles'];tags=a['facet_tags']
    u=a['velocity_m_s'];flow=pv.read(field)
    assert np.array_equal(flow.points,x) and np.array_equal(flow['Velocity'],u)
    assert np.array_equal(np.sort(flow.cells_dict[pv.CellType.TETRA],axis=1),np.sort(tetra,axis=1))
    gradient=core.p1_gradients(x,tetra,u);d,norm,shear=rates(gradient)
    # Independent expanded Cartesian expression for this definition.
    expanded=np.sqrt(2*(gradient[:,0,0]**2+gradient[:,1,1]**2+gradient[:,2,2]**2)
        +(gradient[:,0,1]+gradient[:,1,0])**2+(gradient[:,0,2]+gradient[:,2,0])**2
        +(gradient[:,1,2]+gradient[:,2,1])**2)
    assert np.isfinite(shear).all() and np.allclose(shear,expanded,rtol=5e-15,atol=1e-10)
    angle=.731;q=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    rotated=q@gradient[:1000]@q.T
    rotation_error=float(np.max(abs(rates(rotated)[2]-shear[:1000])))
    assert np.allclose(rates(rotated)[2],shear[:1000],rtol=5e-15,atol=1e-10)
    determinant=np.linalg.det(x[tetra[:,1:]]-x[tetra[:,:1]])
    assert np.all(determinant>0);volume=determinant/6
    owners=core.boundary_owners(tetra,tri)
    center,area,normal=core.wall_geometry(x,tetra,tri,owners)
    surface,ids=core.surface(x,tri)
    flow.cell_data['StrainRateTensor_s_inv']=d.reshape(-1,9)
    flow.cell_data['StrainRateFrobenius_s_inv']=norm
    flow.cell_data['ShearRate_s_inv']=shear
    flow.cell_data['VelocityDivergence_s_inv']=np.trace(gradient,axis1=1,axis2=2)
    flow.cell_data['TetraVolume_m3']=volume
    flow.cell_data['OriginalTetraID_zero_based']=np.arange(len(tetra))
    for name,values in [('StrainRateFrobenius_s_inv',norm),('ShearRate_s_inv',shear)]:
        surface.cell_data[name]=values[owners]
    surface.cell_data['Parent_tetra_zero_based']=owners
    surface.cell_data['Boundary_tag']=tags
    surface.cell_data['Area_m2']=area
    surface.cell_data['Outward_normal']=normal
    flow.save(out/'data/strain_shear_rate_volume_si.vtu')
    surface.save(out/'data/shear_rate_exterior_si.vtp')
    rows=[]
    for name,values in [('ShearRate_s_inv',shear),('StrainRateFrobenius_s_inv',norm)]:
        rows.append(summary(values,volume,'whole_volume',name))
        rows.append(summary(values[owners],area,'whole_exterior',name))
        for tag,label in [(1,'wall'),(4,'inlet_cap'),(3,'O1_cap'),(5,'O2_cap'),(2,'O3_cap')]:
            sel=tags==tag;rows.append(summary(values[owners[sel]],area[sel],label,name))
    with (out/'data/rate_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    # Wall comparison is diagnostic, not the method used to compute shear rate.
    wall=pv.read(ROOT/'input_data/field_diagnostics/data/wall_wss_si.vtp')
    assert np.array_equal(wall.cell_data['Parent_tetra_zero_based'],owners[tags==1])
    mu=compute['material']['mu_Pa_s']
    diff=shear[owners[tags==1]]-wall.cell_data['WSS_raw_Pa']/mu
    maxval=float(shear.max());upper=float(np.ceil(maxval/4000)*4000)
    for name,digest in lock.items():assert sha(ROOT/name)==digest,name
    result=dict(all_pass=True,source_case=compute['case'],source_field_sha256=sha(field),
        source_arrays_sha256=sha(source),gradient_core_sha256=sha(HERE/'audit/wss_core.py'),
        script_sha256=sha(Path(__file__)),source_files_unchanged=True,source_lock=lock,
        definition=dict(G='du_i/dx_j',D='(G+G^T)/2',strain_rate_frobenius='sqrt(D:D)',
            equivalent_shear_rate='sqrt(2*D:D)',unit='s^-1',deviatoric_trace_removal=False),
        material=compute['material'],model='Existing constant-viscosity Newtonian FEM solution; no model change',
        tetrahedra=len(tetra),boundary_triangles=len(tri),wall_triangles=int(sum(tags==1)),
        gradient_method='Exact derivative of P1 tetrahedral velocity; not gradient of speed',
        exterior_mapping='Each exterior triangle receives the scalar from its unique parent tetrahedron',
        surface_display='Raw cell values including inlet/outlet caps; no nodal averaging or scalar interpolation',
        volume_data_saved=True,CFD_calls=0,new_trajectory_integrations=0,
        analytical_checks=checks,expanded_formula_max_error_s_inv=float(np.max(abs(shear-expanded))),
        rotation_invariance_max_error_s_inv=rotation_error,
        wall_shear_rate_minus_WSS_over_mu_max_abs_s_inv=float(np.max(abs(diff))),
        divergence_max_abs_s_inv=float(abs(np.trace(gradient,axis1=1,axis2=2)).max()),
        numerical_divergence_retained=True,mesh_convergence_not_assessed=True,
        scalar_clipping=False,colorbar_range_s_inv=[0,upper],colorbar_ticks_s_inv=np.arange(0,upper+1,4000).tolist(),
        summary=rows,formula_reference='https://cpp.openfoam.org/v13/strainRateViscosityModel_8C_source.html',
        files={str(p.relative_to(out)):sha(p) for p in sorted((out/'data').glob('*'))})
    (out/'COMPUTE_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(all_pass=True,checks=checks,summary=rows[:3],colorbar=[0,upper])),flush=True)

if __name__=='__main__':main()
