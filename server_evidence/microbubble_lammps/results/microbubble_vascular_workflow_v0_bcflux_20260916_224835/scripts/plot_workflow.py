#!/usr/bin/env python3
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
S=Path(__file__).resolve().parents[1];V=S/'visualization';V.mkdir(exist_ok=True)
names=['LOW','MEDIUM','CAPACITY_STRESS_LIMIT','FLUX_WEIGHTED'];colors=['#0072B2','#E69F00','#D55E00','#009E73'];labels=['LOW · source','MEDIUM · admitted','CAPACITY STRESS · admitted','VALID SUPPORT FLUX WEIGHTED']
def read(p):
    with p.open() as f:return list(csv.DictReader(f))
def col(r,key):return np.array([float(x[key]) for x in r])
cases={n:read(S/'runs'/f'{n}_release_mpi1/flux_timeseries.csv') for n in names}
plt.rcParams.update({'font.size':10,'figure.dpi':120,'savefig.dpi':170,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name,title):
    fig.suptitle(title+'\nNOT EXPERIMENTAL CONCENTRATION',fontsize=13);fig.tight_layout(rect=(0,0,1,.92));fig.savefig(V/name);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(12,8))
for ax,n,c,label in zip(axes.flat,names,colors,labels):
    r=cases[n];t=col(r,'time');ax.plot(t,col(r,'target_cumulative'),color='black',label='Cumulative target');ax.step(t,col(r,'admitted_cumulative'),where='post',color=c,label='Admitted');ax.step(t,col(r,'N_source_drawn'),where='post',linestyle=':',color='#777777',label='Source draws');ax.set(title=label,xlabel='Time (s)',ylabel='Bubble count');ax.legend(fontsize=8)
save(fig,'target_vs_actual_injection.png','Cumulative flux control — authoritative Q, independent position support')
for field,file,title in [('active_count','active_bubbles_vs_time.png','Active population — no constant-N feedback'),('pending_count','pending_backlog_vs_time.png','FIFO backlog — fixed size and ID retained')]:
    fig,axes=plt.subplots(2,2,figsize=(12,7))
    for ax,n,c,label in zip(axes.flat,names,colors,labels):
        r=cases[n];ax.step(col(r,'time'),col(r,field),where='post',color=c);ax.set(title=label,xlabel='Time (s)',ylabel='Bubbles');ax.set_ylim(bottom=0)
        if field=='pending_count' and 'STRESS' in n:ax.axhline(20,color='black',ls='--',label='Capacity gate: >20');ax.legend(fontsize=8)
    save(fig,file,title)
fig,axes=plt.subplots(1,3,figsize=(12,4))
for k,ax in enumerate(axes):
    for n,c in zip(names,colors):ax.step(col(cases[n],'time'),col(cases[n],f'outlet_{k}_cumulative'),where='post',color=c,label=n.replace('_LIMIT',''))
    ax.set(title=f'OUTLET_{k}: no natural exit',xlabel='Time (s)',ylabel='Cumulative exits',ylim=(-.1,1));ax.legend(fontsize=7)
save(fig,'outlet_counts_vs_time.png','Actual inlet-transport runs: zero outlet events; cap fixtures are separate tests')
fig,axes=plt.subplots(2,2,figsize=(12,7))
for ax,n,label in zip(axes.flat,names,labels):
    e=read(S/'runs'/f'{n}_release_mpi1/injection_events.csv');draw=[2e6*float(r['radius_m']) for r in e if r['event']=='SOURCE_DRAWN'];adm=[2e6*float(r['radius_m']) for r in e if r['event']=='ADMITTED'];bins=np.arange(.5,5.6,.25);ax.hist(draw,bins,histtype='step',lw=2,label=f'Source (n={len(draw)})');ax.hist(adm,bins,alpha=.5,label=f'Admitted (n={len(adm)})');ax.axvline(2*1.3701249863980327,color='black',ls='--',label='Geometric upper diameter');ax.set(title=label,xlabel='Diameter (µm)',ylabel='Count');ax.legend(fontsize=8)
save(fig,'source_vs_admitted_size_distribution.png','Frozen SonoVue draws and admitted samples — small software-test populations')
fig,axes=plt.subplots(2,2,figsize=(12,7))
for ax,n,label in zip(axes.flat,names,labels):
    r=cases[n][-1];source=float(r['N_source_drawn']);parts=[float(r[k]) for k in ['size_rejected_cumulative','pending_count','active_count']];ax.bar(['Source drawn'],[source],color='#777777');bottom=0
    for value,label2,c in zip(parts,['Size rejected','Pending','Active'],['#CC79A7','#E69F00','#0072B2']):ax.bar(['Resolved states'],[value],bottom=bottom,color=c,label=label2);bottom+=value
    ax.set(title=label,ylabel='Bubbles');ax.legend(fontsize=8);ax.text(.5,.02,'Exits = 0; no deleted or unaccounted bubble',ha='center',transform=ax.transAxes,fontsize=8)
save(fig,'particle_accounting.png','Source = rejected + pending + admitted; admitted = active + exits + terminal')
geom=dict(np.load(S/'geometry/GEOMETRY_ARRAYS.npz'));points=geom['points']*1e6;tri=points[geom['faces']];wallcent=tri[geom['classes']<=2].mean(1);rng=np.random.default_rng(1);wallcent=wallcent[rng.choice(len(wallcent),7000,replace=False)]
fig=plt.figure(figsize=(13,6));ax=fig.add_subplot(121,projection='3d');zoom=fig.add_subplot(122,projection='3d');ax.scatter(*wallcent.T,s=.35,c='#999999',alpha=.15);allx=[]
for n,c in zip(names,colors):
    r=read(S/'runs'/f'{n}_release_mpi1/TRAJECTORIES.csv');ids=sorted(set(int(x['particle_id']) for x in r))
    for j,id in enumerate(ids):
        rr=[x for x in r if int(x['particle_id'])==id];x=np.array([[float(y[k]) for k in ['x_m','y_m','z_m']] for y in rr])*1e6;allx.append(x)
        for a in [ax,zoom]:a.plot(*x.T,color=c,lw=1.3,label=n.replace('_LIMIT','') if j==0 else None)
for a in [ax,zoom]:a.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)');a.view_init(22,-55);a.legend(fontsize=7,loc='upper left')
ax.set_title('Vascular geometry and actual inlet trajectories');zoom.set_title('Trajectory detail; terminal wall stalls retained')
save(fig,'trajectories_3d.png','Frozen-field transport with finite-radius hard wall safety')
section=dict(np.load(S/'geometry/INJECTION_SECTION.npz'));poly=(section['points']-section['center'])@np.c_[section['axis1'],section['axis2']]*1e6
native=read(S/'raw/POSITION_SAMPLES.csv');fig,axes=plt.subplots(1,2,figsize=(12,5))
for ax,mode in zip(axes,['VALID_AREA_UNIFORM','VALID_SUPPORT_FLUX_WEIGHTED']):
    samples=[r for r in native if r['mode']==mode];p=np.array([[float(r[k]) for k in ['x_m','y_m','z_m']] for r in samples]);xy=(p-section['center'])@np.c_[section['axis1'],section['axis2']]*1e6;ax.scatter(*xy.T,s=1,alpha=.15,label='3000 fixed-median-radius position probes')
    for n,c in zip(names,colors):
        if (n=='FLUX_WEIGHTED')!=(mode!='VALID_AREA_UNIFORM'):continue
        e=[r for r in read(S/'runs'/f'{n}_release_mpi1/injection_events.csv') if r['event']=='ADMITTED'];p=np.array([[float(r[k]) for k in ['x_m','y_m','z_m']] for r in e]);xy=(p-section['center'])@np.c_[section['axis1'],section['axis2']]*1e6;ax.scatter(*xy.T,s=35,c=c,edgecolors='black',lw=.5,label=n.replace('_LIMIT','')+' admissions')
    loop=np.vstack([poly,poly[0]]);ax.plot(*loop.T,c='black',lw=1);ax.set_aspect('equal');ax.set(title=mode,xlabel='Plane axis 1 (µm)',ylabel='Plane axis 2 (µm)');ax.legend(fontsize=7)
save(fig,'injection_positions.png','VALID support only + radius clearance — WORKFLOW V0 APPROXIMATION')
# Configurable rolling windows are a postprocessed discrete count diagnostic.
for n in names:
    D=S/'runs'/f'{n}_release_mpi1';r=cases[n];events=read(D/'injection_events.csv');times=np.array([float(e['time_s']) for e in events if e['event']=='ADMITTED']);window=.05
    with (D/'flux_rolling_window.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['time_s','window_s','target_number_flux','actual_admitted_number_flux','outlet_number_flux','relative_flux_error','disclaimer'])
        for x in r[1:]:
            t=float(x['time']);start=max(0,t-window);width=t-start;rate=((times>start)&(times<=t)).sum()/width;target=float(x['target_number_flux']);writer.writerow([t,width,target,rate,0,(rate-target)/target,'NOT EXPERIMENTAL CONCENTRATION'])
print('PLOTS_AND_ROLLING_WINDOWS_CREATED')
