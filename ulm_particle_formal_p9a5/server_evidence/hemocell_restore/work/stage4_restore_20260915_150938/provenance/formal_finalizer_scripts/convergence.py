"""Pure physical-window analysis and confirmation state machine.

Reads observations only; no solver state, population, BC or parameter writes.
"""
import math
import numpy as np
from remote_common import PORTS, OFFSETS, fields, field_path

COLUMNS=['iteration','physical_time','eligible','short_mean_Qin','long_mean_Qin','R_inlet','R_flow','R_CV',
 'R_velocity','R_pressure','flow_fraction_drift','mean_Qout01','mean_Qout02','mean_Qout03',
 'cross_plane_flux_spread','quadrature_error','source_R_velocity','all_gates_pass','consecutive_pass_count',
 'first_converged_iteration','candidate_iteration','confirmation_iteration','failed_gates']

def exact_window(t,y,start,end):
    """Piecewise-linear endpoint interpolation + exact trapezoids on actual times."""
    t=np.asarray(t);y=np.asarray(y)
    if t.ndim!=1 or y.shape[0]!=len(t) or len(t)<2 or np.any(np.diff(t)<=0):raise ValueError('Invalid sample times')
    if start<t[0] or end>t[-1] or end<=start:raise ValueError('Incomplete physical window')
    if not np.isfinite(y).all():raise ValueError('Nonfinite window observations')
    inside=(t>start)&(t<end);tt=np.r_[start,t[inside],end]
    if y.ndim==1:
        yy=np.r_[np.interp(start,t,y),y[inside],np.interp(end,t,y)]
        return float(np.sum(.5*(yy[:-1]+yy[1:])*np.diff(tt))/(end-start))
    return np.array([exact_window(t,y[:,j],start,end) for j in range(y.shape[1])])

def initial_state():
    return {'last_evaluation':None,'first_converged_iteration':None,'candidate_iteration':None,
            'confirmation_iteration':None,'reference_fractions':None,'consecutive_pass_count':0}

def field_residual(run,step,c):
    now=fields(field_path(run,step));old=fields(field_path(run,step-c['formulas']['field_delta_steps']))
    if len(now)!=182694 or not np.array_equal(now['index'],old['index']):raise ValueError('Field node order/domain changed')
    from pathlib import Path
    with np.load(Path(run)/'contracts/control_volume_nodes.npz') as mask:
        if not np.array_equal(np.sort(now['index']),mask['flat_index']):raise ValueError('Field is not frozen physical-fluid domain')
    scale=c['dx_effective_m']/c['dt_s'];punit=c['rho_kg_m3']*scale**2/3
    u=now['u']*scale;uold=old['u']*scale;p=(now['rho']-1)*punit;pold=(old['rho']-1)*punit
    ru=float(np.linalg.norm(u-uold)/max(np.linalg.norm(u),math.sqrt(len(u))*1e-12))
    rp=float(np.linalg.norm(p-pold)/max(np.linalg.norm(p),math.sqrt(len(p))*1.))
    mean=float(np.linalg.norm(u,axis=1).mean());old_mean=float(np.linalg.norm(uold,axis=1).mean())
    return ru,rp,abs(mean-old_mean)/max(abs(mean),1e-12)

def evaluate_point(data,step,c,state,field_values):
    """No IO: window metrics, gates and candidate reset/confirmation rules."""
    state=dict(state);dt=c['dt_s'];end=step*dt;qtarget=c['Qtarget_m3_s'];eps=qtarget*1e-12
    row={k:None for k in COLUMNS};row.update(iteration=step,physical_time=end,eligible=False,
      all_gates_pass=False,consecutive_pass_count=state['consecutive_pass_count'])
    if len(data)==0 or int(data[-1]['iteration'])!=step:raise ValueError('Evaluation must end at a recorded sample')
    if np.any(np.diff(data['iteration'])!=100) or data[0]['iteration']!=0:raise ValueError('Missing scalar samples')
    if not np.allclose(data['time_s'],data['iteration']*dt,rtol=2e-15,atol=1e-20):raise ValueError('Sample time/dt mismatch')
    for name in data.dtype.names:
        if name=='control_volume_mass_balance':
            if len(data)>1 and not np.isfinite(data[name][1:]).all():raise ValueError('Nonfinite CV residual')
        elif not np.isfinite(data[name]).all():raise ValueError('Nonfinite scalar '+name)
    if step<c['execution']['first_possible_evaluation_step'] or end-c['windows']['long_s']<=10*dt:
        state['last_evaluation']=step
        row['failed_gates']='NOT_EVALUATED_INCOMPLETE_PHYSICAL_WINDOWS'
        return row,state,{'windows':{},'gates':{},'stop':False,'new_candidate':False}
    if field_values is None:raise ValueError('Eligible evaluation missing independent field residual')
    row['eligible']=True;windows={};rhoscale=c['rho_kg_m3'];t=data['time_s']
    flowcols=[f'{p}_{k}dx' for p in PORTS for k in OFFSETS]
    volume=np.column_stack([data[n] for n in flowcols])
    mass=np.column_stack([data[f'mass_outward_g{g}'] for g in [2,8,14,20]])
    for name in ['short','long']:
        duration=c['windows'][name+'_s'];start=end-duration
        means=exact_window(t,volume,start,end).reshape(4,3);primary=means[:,1]
        mean_mass=exact_window(t,mass,start,end)
        mass_start=float(np.interp(start,t,data['control_volume_mass']))
        mass_end=float(data['control_volume_mass'][-1]);derivative=(mass_end-mass_start)/duration
        denominator=max(abs(mean_mass[0]),float(np.abs(mean_mass[1:]).sum()),rhoscale*eps)
        rcv=abs(derivative+float(mean_mass.sum()))/denominator
        mean_abs_q=exact_window(t,np.abs(data['Qin_4dx']),start,end)
        rf=abs(primary[0]-float(primary[1:].sum()))/max(mean_abs_q,eps)
        ri=abs(primary[0]-qtarget)/qtarget
        total=float(primary[1:].sum());fractions=(primary[1:]/total).tolist() if total>0 else None
        spread=(np.max(means,axis=1)-np.min(means,axis=1))/np.maximum(np.mean(np.abs(means),axis=1),eps)
        windows[name]={'Qin':float(primary[0]),'Qout':primary[1:].tolist(),'R_flow':rf,'R_inlet':ri,'R_CV':rcv,
          'source_R_mass':abs(primary[0]-total)/max(abs(primary[0]),np.finfo(float).tiny),
          'fractions':fractions,'all_plane_means':means.tolist(),'cross_plane_spread':spread.tolist(),
          'control_mass_start':mass_start,'control_mass_end':mass_end,'dM_dt':derivative,'mean_mass_outward':mean_mass.tolist(),
          'mean_speed_m_s':exact_window(t,data['mean_speed_m_s'],start,end),
          'significant_backflow':[bool(v<0 and abs(v)>.05*abs(primary[0])) for v in primary[1:]],
          'start_s':start,'end_s':end,'mean_method':'exact physical endpoints, piecewise linear trapezoidal mean'}
    fs,fl=windows['short']['fractions'],windows['long']['fractions'];drift=None
    if fs is not None and fl is not None:
        drift=float(np.max(np.abs(np.array(fs)-fl)))
        if state['reference_fractions'] is not None:
            for name in ['short','long']:
                drift=max(drift,float(np.max(np.abs(np.array(windows[name]['fractions'])-state['reference_fractions'][name]))))
    quad_errors=[]
    for p in PORTS:
        for k in OFFSETS:
            coarse=float(data[f'{p}_{k}dx'][-1]);fine=float(data[f'{p}_{k}dx_dx2'][-1])
            quad_errors.append(abs(coarse-fine)/max(abs(fine),eps))
    ru,rp,source_ru=field_values
    row.update(short_mean_Qin=windows['short']['Qin'],long_mean_Qin=windows['long']['Qin'],
      R_inlet=max(w['R_inlet'] for w in windows.values()),R_flow=max(w['R_flow'] for w in windows.values()),
      R_CV=max(w['R_CV'] for w in windows.values()),R_velocity=ru,R_pressure=rp,
      flow_fraction_drift=drift,cross_plane_flux_spread=max(max(w['cross_plane_spread']) for w in windows.values()),
      quadrature_error=max(quad_errors),source_R_velocity=source_ru)
    for i in range(3):row[f'mean_Qout0{i+1}']=windows['long']['Qout'][i]
    mapping={'R_inlet':'R_inlet','R_flow':'R_flow','R_CV':'R_CV','R_velocity':'R_velocity','R_pressure':'R_pressure',
      'flow_fraction_drift':'flow_fraction_drift','cross_plane_flux_spread':'cross_plane_spread','quadrature_error':'quadrature_dx_vs_dx2'}
    gates={key:bool(row[key] is not None and math.isfinite(row[key]) and row[key]<=c['gates'][gate]) for key,gate in mapping.items()}
    for i in range(3):gates[f'outlet_0{i+1}_mean_direction']=all(w['Qout'][i]>0 for w in windows.values())
    gates['positive_mean_inlet']=all(w['Qin']>0 for w in windows.values())
    passed=all(gates.values());row['all_gates_pass']=passed;new_candidate=False;stop=False
    if not passed:
        state.update(candidate_iteration=None,reference_fractions=None,consecutive_pass_count=0)
    elif state['candidate_iteration'] is None:
        state['candidate_iteration']=step
        if state['first_converged_iteration'] is None:state['first_converged_iteration']=step
        state['reference_fractions']={name:windows[name]['fractions'] for name in windows}
        state['consecutive_pass_count']=1;new_candidate=True
    elif step-state['candidate_iteration']>=c['execution']['confirmation_separation_steps']:
        state['confirmation_iteration']=step;state['consecutive_pass_count']=2;stop=True
    state['last_evaluation']=step
    row.update(consecutive_pass_count=state['consecutive_pass_count'],first_converged_iteration=state['first_converged_iteration'],
      candidate_iteration=state['candidate_iteration'],confirmation_iteration=state['confirmation_iteration'],
      failed_gates=';'.join(name for name,passed in gates.items() if not passed))
    # Source companion range across instantaneous a,b,c, kept distinct from mean drift.
    companion=[]
    for when in [end-c['windows']['long_s'],end-c['windows']['short_s'],end]:
        oq=np.array([np.interp(when,t,data[f'Qout0{i}_4dx']) for i in [1,2,3]])
        companion.append((oq/oq.sum()).tolist() if oq.sum()>0 else None)
    source_drift=None if any(x is None for x in companion) else float(np.ptp(companion,axis=0).max())
    return row,state,{'windows':windows,'gates':gates,'quadrature_errors':quad_errors,
                     'source_fraction_drift_companion':source_drift,'stop':stop,'new_candidate':new_candidate}

def check_observed_safety(data,c):
    bounds=c['runtime_safety']['fluid_density_bounds'];mach=c['runtime_safety']['mach_max']
    if np.any(data['rho_min']<bounds[0]) or np.any(data['rho_max']>bounds[1]):raise ValueError('Observed density safety failure')
    if np.any(data['Mach_max']>mach) or np.any(data['ghost_Mach_max']>mach) or np.any(data['ghost_rho_min']<=0):raise ValueError('Observed Mach/ghost safety failure')
