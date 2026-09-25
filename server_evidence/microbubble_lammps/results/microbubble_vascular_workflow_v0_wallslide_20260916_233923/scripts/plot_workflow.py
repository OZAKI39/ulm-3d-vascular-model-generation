from pathlib import Path
import csv,json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
S=Path(__file__).resolve().parents[1];B=S.parent/'bcflux_20260916_224835';O=S/'visualization';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.grid':True,'figure.dpi':140})
def rows(p):return [r for r in csv.DictReader(p.open()) if r.get('disclaimer')=='NOT EXPERIMENTAL CONCENTRATION' and None not in r and all(v is not None for v in r.values())]
def arr(r,keys):return np.array([[float(x[k]) for k in keys] for x in r])
def save(fig,name,title):
 fig.suptitle(title);fig.text(.5,.015,'TEST_ONLY FLUX — NOT EXPERIMENTAL CONCENTRATION\nGEOMETRIC CONSTRAINT ONLY | ENGINEERING_TRANSIENT_FIELD_ONLY',ha='center',fontsize=8)
 fig.tight_layout(rect=[0,.075,1,.94]);fig.savefig(O/name);plt.close(fig)
def run(case):return S/'runs'/f'{case}_release_mpi1'
def trajectory(ax,r,label='',style='-'):
 for id in sorted(set(x['particle_id'] for x in r),key=int):
  q=[x for x in r if x['particle_id']==id];q=q[::max(1,len(q)//5000)]+q[-1:];xyz=arr(q,['x_m','y_m','z_m'])*1e6
  ax.plot(*xyz.T,style,lw=1.2,label=label+' ID '+id)
 ax.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)');ax.legend(fontsize=7)
for c in ['LOW','MEDIUM']:
 fig=plt.figure(figsize=(9,6));ax=fig.add_subplot(projection='3d');trajectory(ax,rows(B/'runs'/f'{c}_release_mpi1/TRAJECTORIES.csv'),'Old','--');trajectory(ax,rows(run(c)/'TRAJECTORIES.csv'),'New')
 save(fig,f'old_vs_new_{c}_trajectory.png',c+': original stall replay and constrained trajectory')
wall={c:rows(run(c)/'wall_constraint_timeseries.csv') for c in ['LOW','MEDIUM']}
fig,axs=plt.subplots(2,1,figsize=(9,6),sharex=False)
for ax,c in zip(axs,wall):
 for id in sorted(set(r['particle_id'] for r in wall[c]),key=int):
  q=[r for r in wall[c] if r['particle_id']==id and r['stage']=='1'];v=arr(q,['time','gap']);ax.plot(v[:,0],v[:,1]*1e9,label='ID '+id)
 ax.axhline(.1,color='k',ls=':',label='safety margin');ax.set(title=c,ylabel='surface gap (nm)',yscale='log',xlabel='time (s)');ax.legend(fontsize=7,ncol=4)
save(fig,'gap_vs_time_wall_contact.png','Wall gap at hydrodynamic evaluation points')
fig,axs=plt.subplots(2,1,figsize=(9,6))
for ax,c in zip(axs,wall):
 q=[r for r in wall[c] if r['wall_constraint_active']=='1'];v=arr(q,['time','vn_raw','vn_used']);ax.scatter(v[:,0],v[:,1]*1e6,s=5,label='raw');ax.scatter(v[:,0],v[:,2]*1e6,s=5,label='used');ax.set(title=c,xlabel='time (s)',ylabel='normal velocity (µm/s)');ax.legend()
save(fig,'raw_vs_constrained_normal_velocity.png','Only unsafe inward normal velocity is removed')
fig,ax=plt.subplots(figsize=(7,6))
for c in wall:
 q=[r for r in wall[c] if r['wall_constraint_active']=='1'];v=arr(q,['tangential_speed_raw','tangential_speed_used'])*1e6;ax.scatter(v[:,0],v[:,1],s=8,label=c)
lims=ax.get_xlim();ax.plot(lims,lims,'k--',label='unchanged');ax.set(xlabel='raw tangential speed (µm/s)',ylabel='used tangential speed (µm/s)');ax.legend()
save(fig,'tangential_speed_before_after.png','Tangential preservation at every active event')
fig=plt.figure(figsize=(9,6));ax=fig.add_subplot(projection='3d')
for c in wall:
 q=[r for r in wall[c] if r['wall_constraint_active']=='1'];xyz=arr(q,['x_eval','y_eval','z_eval'])*1e6;ax.scatter(*xyz.T,s=6,label=c)
ax.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)');ax.legend();save(fig,'wall_constraint_events_3d.png','Locations of accepted kinematic wall constraints')
fig,axs=plt.subplots(2,1,figsize=(9,6))
for ax,c in zip(axs,['LOW','MEDIUM']):
 for root,label in [(B,'old'),(S,'new')]:
  v=arr(rows(root/'runs'/f'{c}_release_mpi1/SOLVER_HISTORY.csv'),['time_s','dt_s']);ax.semilogy(v[:,0],v[:,1],label=label)
 ax.set(title=c,xlabel='time (s)',ylabel='accepted dt (s)');ax.legend()
save(fig,'dt_history_old_vs_new.png','Accepted timestep histories; terminal reason retained separately')
long=rows(run('LONG_TRANSPORT')/'TRAJECTORIES.csv');evidence=json.loads((S/'validation/LONG_TRANSPORT_ABORT_EVIDENCE.json').read_text())['cases'][0];state={'time_s':evidence['last_flushed_accounting_time_s'],'reason':'STOP_WALL_NORMAL_AMBIGUITY','outlet_exits':evidence['outlet_exits']}
fig=plt.figure(figsize=(9,6));ax=fig.add_subplot(projection='3d');trajectory(ax,long);save(fig,'LONG_TRANSPORT_trajectories.png',f"LONG_TRANSPORT: t={state['time_s']:.6f} s; {state['reason']}")
fig,ax=plt.subplots(figsize=(9,5))
for c in ['LOW','MEDIUM','LONG_TRANSPORT']:
 v=arr(rows(run(c)/'flux_timeseries.csv'),['time','active_count']);ax.step(v[:,0],v[:,1],where='post',label=c)
ax.set(xlabel='time (s)',ylabel='active bubbles');ax.legend();save(fig,'active_bubbles_vs_time.png','Adaptive injection remains enabled; actual particle accounting')
fig,ax=plt.subplots(figsize=(9,5));flux=rows(run('LONG_TRANSPORT')/'flux_timeseries.csv');v=arr(flux,['time','outlet_0_cumulative','outlet_1_cumulative','outlet_2_cumulative'])
for i in range(3):ax.step(v[:,0],v[:,i+1],where='post',label=f'OUTLET_{i}')
ax.set(xlabel='time (s)',ylabel='natural outlet exits');ax.legend()
if state['outlet_exits']==0:ax.text(.5,.5,'NATURAL_OUTLET_STATUS = NOT_OBSERVED',transform=ax.transAxes,ha='center',bbox=dict(facecolor='white',alpha=.9))
save(fig,'outlet_events_long_transport.png','Natural events only; synthetic cap fixtures excluded')
fig,ax=plt.subplots(figsize=(10,6));ax.set_aspect('equal');ax.plot([-2,3],[0,0],color='black',lw=3);ax.fill_between([-2,3],-.5,0,color='lightgray');circle=plt.Circle((0,1),.32,fc='white',ec='navy',lw=2);ax.add_patch(circle)
for end,color,label,pos in [((1.5,.3),'tab:orange','V_raw',(1.25,.45)),((1.5,1),'tab:blue','V_used = V_tangent',(1.35,1.15)),((0,.3),'tab:red','V_inward removed',(-1.1,.4)),((0,1.8),'tab:green','n: wall → center',(-.9,1.85))]:
 ax.annotate('',xy=end,xytext=(0,1),arrowprops=dict(arrowstyle='->',color=color,lw=2));ax.text(*pos,label,color=color,ha='center',fontsize=11)
ax.text(.4,2.4,'V_raw = V_tangent + V_inward\nV_used = V_raw − min(V_raw · n, 0)n\nOmega_used = Omega_raw',ha='center',fontsize=12)
ax.text(.5,-.8,'GEOMETRIC CONSTRAINT ONLY\nNO WALL FORCE\nNO WALL HYDRODYNAMICS',ha='center',weight='bold',fontsize=12)
ax.set(xlim=(-2,3),ylim=(-1.3,3.1));ax.axis('off');save(fig,'kinematic_wall_sliding_scheme.png','Activated only when the raw stage segment is unsafe')
syn=json.loads((S/'validation/SYNTHETIC_WALL_VALIDATION.json').read_text());fig,axs=plt.subplots(1,2,figsize=(10,4.5))
for c in syn['curved']:
 r=c['refinement'];dt=[x['dt_max'] for x in r];axs[0].loglog(dt,[x['trajectory_error_m'] for x in r],'o-',label=c['shape']);axs[1].loglog(dt,[x['max_correction_m'] for x in r],'o-',label=c['shape'])
for ax,title in zip(axs,['trajectory error (m)','maximum correction (m)']):ax.set(xlabel='dt_max (s)',ylabel=title);ax.legend()
save(fig,'synthetic_curved_refinement.png','Cylinder and sphere: refinement of trajectories and corrections')
print('PLOTS_WRITTEN',len(list(O.glob('*.png'))))
