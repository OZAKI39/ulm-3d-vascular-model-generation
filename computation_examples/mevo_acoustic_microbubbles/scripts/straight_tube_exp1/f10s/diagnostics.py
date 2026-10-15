"""Finite force-stability and Y-scope diagnostics, no particle trajectories."""
import copy
import time
from dataclasses import replace
import numpy as np
from scipy.optimize import minimize_scalar, brentq
from scipy.linalg import expm
from mevo_mb.acoustics.water_disc_proxy import WaterDiscs
from mevo_mb.fast_sim.nearfield import pair_reference
from mevo_mb.fast_sim.f10g_contact import resistance
from .common import *

def field(c,x,u,discs=None):
    pb,gb=(discs or c['discs']).evaluate(x)
    return u@pb,np.einsum('i,inj->nj',u,gb)

def primary(c,x,u,discs=None):
    pp,gg=field(c,x,u,discs)
    v=c['response'].query(pp)
    return -.5*np.real(v[:,0,None]*gg.conj()),pp,gg,v

def jacobian(fn,x,h):
    return np.column_stack([(fn(x+np.eye(3)[i]*h)-fn(x-np.eye(3)[i]*h))/(2*h) for i in range(3)])

def frozen(c,x,u):
    f,p,g,v=primary(c,x,u)
    phys=c['physics'];chi=c['response'].query(np.array([complex(c['response'].knots[0])]))[0,0]/c['response'].knots[0]
    secondary=pair_reference(x,v,p,phys,chi)
    far=pair_reference(x,v,p,replace(phys,force_model='f1_smooth_near_field'),chi)
    op=resistance(x,phys.radius_m);zeta=6*np.pi*phys.viscosity_Pa_s*phys.radius_m
    T=np.tile(np.eye(3),(len(x),1))
    rigid_v=np.linalg.solve(zeta*T.T@op['A']@T,T.T@(f+secondary).ravel())
    free_v=np.linalg.solve(zeta*op['A'],(f+secondary).ravel()).reshape(-1,3)
    i,j=np.triu_indices(len(x),1);delta=x[j]-x[i];dist=np.linalg.norm(delta,axis=1)
    gap_rate=np.einsum('ij,ij->i',free_v[j]-free_v[i],delta/dist[:,None])
    edges=dist<=2.1*phys.radius_m
    return dict(positions_m=x,pressure_Pa=p,gradient_Pa_m=g,primary_N=f,secondary_N=secondary,F2_increment_N=secondary-far,
        COM_velocity_m_s=rigid_v,unconstrained_particle_velocity_m_s=free_v,edge_pairs=np.column_stack([i[edges],j[edges]]),
        edge_gap_rates_m_s=gap_rate[edges],max_free_approach_1us_m=float(max(0,-gap_rate[edges].min())*1e-6),
        compression_virial_N_m=float(np.sum((x-x.mean(0))*(f+secondary))),
        scope='FROZEN_SHAPE_COM_DIAGNOSTIC; T^T zeta A T U = T^T F. Pair lubrication cancels for uniform translation; internal pair forces sum to zero. No contact-step admission inferred.')

def traps(c,cfg):
    start=time.perf_counter(); records=[];curves={};x0=c['x0'];h=cfg['jacobian_delta_m']
    for amp in cfg['amplitudes_Pa_equivalent']:
      for deg in cfg['phases_deg']:
        if time.perf_counter()-start>100:raise RuntimeError('C budget exhausted')
        u=np.array([0,0,amp,amp*np.exp(1j*np.deg2rad(deg))])
        peak=minimize_scalar(lambda x:-abs(field(c,[[x,0,0]],u)[0][0])**2,bounds=(-cfg['peak_search_halfwidth_m'],cfg['peak_search_halfwidth_m']),method='bounded',options={'xatol':1e-14}).x
        row=dict(amplitude_Pa_equivalent=amp,phase_deg=deg,pressure_peak_x_m=peak)
        grid=peak+np.linspace(-cfg['force_search_halfwidth_m'],cfg['force_search_halfwidth_m'],161)
        pts=np.column_stack([grid,np.zeros((len(grid),2))]); pp,gg=field(c,pts,u)
        curve=dict(x_m=grid,pressure_Pa=pp,gradient_Pa_m=gg,phase_rad=np.angle(pp))
        try:
            f,*_=primary(c,pts,u);curve['primary_force_N']=f
            root=brentq(lambda x:primary(c,[[x,0,0]],u)[0][0,0],grid[0],grid[-1],xtol=1e-16)
            point=np.array([root,0,0]);J=jacobian(lambda q:primary(c,[q],u)[0][0],point,h)
            J2=jacobian(lambda q:primary(c,[q],u)[0][0],point,h/2)
            eig=np.linalg.eigvals(J);axial=J[0,0]<0;all3=bool(np.all(eig.real<0))
            row.update(status='PASS',zero_force_x_m=root,peak_force_zero_gap_m=root-peak,force_J_N_m=J,
              force_J_halfstep_max_difference_N_m=float(abs(J-J2).max()),force_J_eigenvalues_N_m=eig.real,
              equilibrium_force_N=primary(c,[point],u)[0][0],axial_restoring=bool(axial),three_dimensional_restoring=all3,
              stability_label='LOCAL_3D_RESTORING' if all3 else ('AXIAL_RESTORING_ONLY' if axial else 'NOT_A_RESTORING_TRAP'))
            row['original_geometry_at_original_COM']=frozen(c,x0,u)
            row['original_geometry_at_single_bubble_equilibrium']=frozen(c,x0+point,u)
            if axial:
                # Fixed internal geometry, nonlinear actual APECSS forces; find group equilibrium separately.
                def velocity(q):
                    f,*_=primary(c,x0+q,u)
                    op=resistance(x0,c['physics'].radius_m)['A'];T=np.tile(np.eye(3),(32,1))
                    zeta=6*np.pi*c['physics'].viscosity_Pa_s*c['physics'].radius_m
                    return np.linalg.solve(zeta*T.T@op@T,T.T@f.ravel())
                center=brentq(lambda a:velocity(np.array([a,0,0]))[0],grid[0],grid[-1],xtol=1e-16)
                q=np.array([center,0,0]);JV=jacobian(velocity,q,h);JV2=jacobian(velocity,q,h/2)
                tau=-1/JV[0,0] if JV[0,0]<0 else None
                shift_per_deg=-.842628e-6
                offset=np.array([-shift_per_deg,0,0]);error=expm(JV*.002)@offset
                row['frozen_COM']=dict(model='FROZEN_SHAPE_COM_DIAGNOSTIC',equilibrium_m=q,
                    drift_J_s_inverse=JV,J_halfstep_max_difference=float(abs(JV-JV2).max()),
                    tau_x_s=tau,eigenvalues_s_inverse=np.linalg.eigvals(JV).real,
                    initial_offset_1um_drift_m_s=velocity(q+np.array([1e-6,0,0])),
                    imposed_peak_step_m=shift_per_deg,predicted_following_displacement_2ms_m=error-offset,
                    prediction='Linearized response to a 1-degree step using F10R measured local peak slope; not a particle trajectory or a phase-ramp measurement',
                    full_moving_trap_run='NOT_RUN',reason='Relaxation much longer than 2ms' if tau and tau>.002 else 'Only instantaneous diagnostic; no full admission')
        except ValueError as ex:
            row.update(status='REJECTED',reason=str(ex),force_stability='NOT_RUN')
        curves[f'A{amp:g}_phase{deg:+g}']=curve;records.append(row)
    write(BASE/'diagnostics/TRAPS.json',dict(records=records,curves=curves,wall_s=time.perf_counter()-start,identity=identity(),particle_trajectories_run=0))
    return records

def yscope(c,cfg):
    start=time.perf_counter();geo=read(OLD/'CONFIG_SNAPSHOT.json')['geometry'];rows=[]
    settings=[('source_distance',a,0,2.5e6,0) for a in cfg['Y_distance_offsets_m']]
    settings += [('COM_y',0,a,2.5e6,0) for a in cfg['COM_y_offsets_m']]
    settings += [('frequency',0,0,a,0) for a in cfg['frequencies_Hz']]
    settings += [('phase',0,0,2.5e6,a) for a in cfg['single_phase_deg']]
    for kind,offset,cy,freq,phase in settings:
        elements=copy.deepcopy(geo['pzt']);elements[0]['center_local_m'][1]-=offset
        d=WaterDiscs(elements,frequency_Hz=freq);u=np.array([cfg['source_Pa_equivalent']*np.exp(1j*np.deg2rad(phase)),0,0,0])
        point=np.array([[0,cy,0.]]);p,g=field(c,point,u,d)
        px,_=field(c,point,np.array([0,0,cfg['source_Pa_equivalent'],0]),d)
        wavelength=1500/freq;zero=(.0035**2-(2*wavelength)**2)/(4*wavelength)
        row=dict(kind=kind,source_distance_offset_m=offset,COM_y_m=cy,frequency_Hz=freq,phase_deg=phase,
            pressure_Pa=p[0],pressure_magnitude_Pa=abs(p[0]),gradient_Pa_m=g[0],Y_over_X_pressure=float(abs(p[0]/px[0])),
            axial_second_zero_distance_m=zero,zero_world_y_m=elements[0]['center_local_m'][1]+zero,
            cavity_legal=bool(c['cavity'].query(point[0],c['physics'].radius_m)['valid']))
        if freq==2.5e6:
            try: row.update(primary_force_N=primary(c,point,u,d)[0][0],force_status='PASS')
            except ValueError as ex:row.update(force_status='REJECTED',reason=str(ex))
        else:row.update(force_status='NOT_RUN_FREQUENCY_OUTSIDE_APECSS',primary_force_N=None)
        rows.append(row)
    write(BASE/'diagnostics/Y_SCOPE.json',dict(rows=rows,wall_s=time.perf_counter()-start,identity=identity(),scope=cfg['scope'],particle_trajectories_run=0))

def main():
    start=time.perf_counter();(BASE/'diagnostics').mkdir(exist_ok=True)
    cfg=read(CONFIG);c=context()
    records=traps(c,cfg['C']);yscope(c,cfg['D'])
    write(BASE/'diagnostics/EXECUTION.json',dict(wall_s=time.perf_counter()-start,identity=identity(),particle_trajectories_run=0,
         frequency_Hz=2.5e6,response_pressure_knots_Pa=c['response'].knots,training_ready=False))
    print([(r['amplitude_Pa_equivalent'],r['phase_deg'],r['status'],r.get('stability_label'),r.get('frozen_COM',{}).get('tau_x_s')) for r in records])
if __name__=='__main__':main()
