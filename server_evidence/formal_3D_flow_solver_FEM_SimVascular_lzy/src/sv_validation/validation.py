"""Fail-closed validation contracts for SV0; pure array tests are not real-run evidence."""
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

def real_boundary_conditions(bcs):
    require(bcs.get('wall')=={'type':'Dirichlet','value':[0.,0.,0.]},'Rigid no-slip wall required')
    require(set(bcs.get('outlets',{}))=={'OUTLET_01','OUTLET_02','OUTLET_03'},'Missing outlet')
    require(all(v=={'type':'Neumann','traction':0.} for v in bcs['outlets'].values()),'All outlets must be zero traction')
    inlet=bcs['inlet'];require(inlet['type']=='Dirichlet' and inlet['profile']=='plug','Primary inlet must be declared plug profile')
    require(inlet.get('equivalent_to_integral_only') is False,'Prescribed profile is not integral-only boundary data')
    return True

def surface_geometry(distances_a,distances_b,h_wall,volume_before,volume_after):
    distances=np.r_[distances_a,distances_b]
    require(len(distances)>0 and np.isfinite(distances).all() and np.all(distances>=0),'Invalid surface distance data')
    require(h_wall>0 and volume_before>0 and volume_after>0,'Positive physical scales required')
    stats={'maximum':float(np.max(distances)),'P95':float(np.quantile(distances,.95)),
           'RMS':float(np.sqrt(np.mean(distances**2))),'relative_volume_error':abs(volume_after-volume_before)/volume_before}
    require(stats['maximum']<=.1*h_wall and stats['P95']<=.03*h_wall and stats['RMS']<=.02*h_wall and stats['relative_volume_error']<=.005,'Surface geometry drift exceeds frozen user gates')
    return stats

def port_geometry(area_relative_error,centroid_shift,equivalent_radius,normal_dot,plane_deviation,planarity_tolerance):
    require(np.isfinite([area_relative_error,centroid_shift,equivalent_radius,normal_dot,plane_deviation,planarity_tolerance]).all(),'Non-finite port geometry')
    require(equivalent_radius>0 and planarity_tolerance>0,'Positive port scales required')
    require(area_relative_error<=.005 and centroid_shift<=.02*equivalent_radius and normal_dot>=.999 and plane_deviation<=planarity_tolerance,'Port geometry gate failed')
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
    require(error<=1e-12,'Written inlet field is not normalized to target Q')
    return error

def fields_finite(velocity,pressure):
    require(np.size(velocity)>0 and np.size(pressure)>0,'Empty solution fields')
    require(np.isfinite(velocity).all() and np.isfinite(pressure).all(),'Non-finite solution fields')
    return True

def result_exists(returncode,path):
    require(returncode==0,'Solver exit code is nonzero')
    p=Path(path);require(p.is_file() and p.suffix=='.vtu' and p.stat().st_size>0,'No actual result VTU')
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
    error=max(abs(solver_flux[n]-independent_flux[n])/target for n in solver_flux)
    require(np.isfinite(error) and error<=1e-6,'POSTPROCESS_DISAGREEMENT')
    return error

def steady_state(velocity_errors,flow_errors):
    eu=np.asarray(velocity_errors);eq=np.asarray(flow_errors)
    require(len(eu)==len(eq) and len(eu)>=5,'At least five saved intervals, hence at least six states, required')
    require(np.isfinite(eu).all() and np.isfinite(eq).all(),'Non-finite steady diagnostics')
    return bool(np.all(eu[-5:]<=1e-5) and np.all(eq[-5:]<=1e-6))
