"""Validation contracts for native mesh and flow data."""
from pathlib import Path
import numpy as np

ROLES={'WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03'}
SI_UNITS={'length':'m','mass':'kg','time':'s','velocity':'m/s','pressure':'Pa','density':'kg/m^3','dynamic_viscosity':'Pa*s','volume_flow':'m^3/s'}

class ValidationError(ValueError):pass

def require(condition,message):
    if not condition:raise ValidationError(message)

def reference_condition(config):
    require(config['condition_type']=='REFERENCE_NUMERICAL_CONDITION' and config['experimental'] is False,'Only frozen numerical reference permitted')
    require(config['warning']=='NOT_EXPERIMENTAL_PUMP_FLOW','Reference provenance warning required')
    require(config['units']==SI_UNITS,'Mixed units or unvalidated scale adapter')
    p=config['physics'];rho=p['density_kg_m3'];nu=p['kinematic_viscosity_m2_s'];mu=p['dynamic_viscosity_pa_s'];Q=p['inlet_volume_flow_m3_s']
    require(np.isfinite([rho,nu,mu,Q]).all() and min(rho,nu,mu,Q)>0,'Positive finite reference properties required')
    require(np.isclose(mu,rho*nu,rtol=1e-14,atol=0),'Dynamic viscosity must be rho*nu')
    return p

def face_mapping(mapping,original_roles):
    require(len(mapping)==5,'Five named exterior roles required')
    require(set(r['role'] for r in mapping)==ROLES,'Exactly one inlet, three outlets, one wall required')
    require(len({r['sv_face_id'] for r in mapping})==5,'No split or merged semantic faces')
    for r in mapping:
        require(original_roles.get(str(r['original_entity_id']))==r['role'],'Role swapped relative to original labels')
    return True

def triangle_flux(points,triangles,velocity):
    p=np.asarray(points,float);t=np.asarray(triangles,int);u=np.asarray(velocity,float)
    require(np.isfinite(p).all() and np.isfinite(u).all(),'Non-finite boundary fields')
    area_vector=.5*np.cross(p[t[:,1]]-p[t[:,0]],p[t[:,2]]-p[t[:,0]])
    # Exact integral of piecewise linear nodal velocity on oriented triangles.
    return float(np.sum(np.einsum('ij,ij->i',u[t].mean(axis=1),area_vector)))

def normalized_inlet(actual_inflow,target):
    require(target>0 and np.isfinite([actual_inflow,target]).all(),'Invalid inlet flow')
    error=abs(actual_inflow-target)/target
    require(error<=1e-10,'Written inlet field is not normalized to target Q')
    return error

def fields_finite(velocity,pressure):
    require(np.size(velocity)>0 and np.size(pressure)>0,'Empty solution fields')
    require(np.isfinite(velocity).all() and np.isfinite(pressure).all(),'Non-finite solution fields')
    return True

def mass_balance(target,inlet_signed,outlets):
    require(set(outlets)=={'OUTLET_01','OUTLET_02','OUTLET_03'},'Need all three measured outlet flows')
    require(target>0 and np.isfinite([target,inlet_signed,*outlets.values()]).all(),'Non-finite flux')
    incoming=-inlet_signed;total=sum(outlets.values());error=abs(incoming-target)/target;closure=abs(total-incoming)/target
    require(error<=1e-6 and closure<=1e-6,'Actual solution fails inlet or mass conservation gate')
    require(total!=0,'Outlet fractions undefined')
    return {'epsilon_Q':error,'epsilon_mass':closure,'Q_in':incoming,'Q_out_total':total,
            'outlet_fractions':{n:q/total for n,q in outlets.items()}}

def integral_agreement(solver_flux,independent_flux,target):
    require(target>0 and set(solver_flux)==set(independent_flux)=={'INLET','OUTLET_01','OUTLET_02','OUTLET_03'},'Boundary integral identities incomplete')
    require(np.isfinite([target,*solver_flux.values(),*independent_flux.values()]).all(),'POSTPROCESS_DISAGREEMENT: non-finite boundary integral')
    error=max(abs(solver_flux[n]-independent_flux[n])/target for n in solver_flux)
    require(np.isfinite(error) and error<=1e-6,'POSTPROCESS_DISAGREEMENT')
    return error

def steady_state(velocity_errors,flow_errors):
    eu=np.asarray(velocity_errors);eq=np.asarray(flow_errors)
    require(len(eu)==len(eq) and len(eu)>=5,'At least five saved intervals, hence at least six states, required')
    require(np.isfinite(eu).all() and np.isfinite(eq).all(),'Non-finite steady diagnostics')
    require(np.all(eu>=0) and np.all(eq>=0),'Steady norm errors must be nonnegative')
    return bool(np.all(eu[-5:]<=1e-5) and np.all(eq[-5:]<=1e-6))


def parse_result_vtu(path,velocity_name='Velocity',pressure_name='Pressure'):
    import pyvista as pv
    mesh=pv.read(path)
    require(mesh.n_cells>0 and mesh.n_points>0,'Empty result VTU')
    require(velocity_name in mesh.point_data and pressure_name in mesh.point_data,'Missing velocity or pressure array')
    u=np.asarray(mesh.point_data[velocity_name]);p=np.asarray(mesh.point_data[pressure_name])
    require(u.shape==(mesh.n_points,3) and p.size==mesh.n_points,'Unexpected result array dimensions')
    fields_finite(u,p)
    return mesh,u,p
