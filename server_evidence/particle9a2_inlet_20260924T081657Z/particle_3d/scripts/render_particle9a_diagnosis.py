#!/usr/bin/env python3
"""Six static review figures, generated only from saved diagnostic evidence."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from analyze_particle9a_diagnosis import DATA,REPORT,load,dump,csvwrite

FIG=REPORT/'figures'
COLORS={'COMPLETED':'#258477','INLET_ESCAPE':'#c33c49','SUBDIVISION_STOP':'#d3871b','TIME_LIMIT':'#7358a5'}
LABELS={'COMPLETED':'Completed','INLET_ESCAPE':'Inlet escape','SUBDIVISION_STOP':'Subdivision stop','TIME_LIMIT':'Time limit'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'axes.labelsize':10,'axes.titlesize':12,'savefig.facecolor':'white','figure.facecolor':'white',
    'axes.facecolor':'white','pdf.fonttype':42,'svg.fonttype':'none'})
S=json.loads((DATA/'analysis_summary.json').read_text())
MANIFEST=[]

def save(fig,name,sources):
    FIG.mkdir(exist_ok=True);fig.savefig(FIG/(name+'.png'),dpi=230,bbox_inches='tight')
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    MANIFEST.append(dict(figure='figures/'+name+'.png',pdf='figures/'+name+'.pdf',sources=sources))

def table1():
    rows=list(csv.DictReader((DATA/'ab_outcomes.csv').open()));ids=sorted({int(r['particle_id']) for r in rows})
    fig,ax=plt.subplots(figsize=(9,7));ax.axis('off')
    body=[];colors=[]
    for pid in ids:
        pair=[next(r for r in rows if int(r['particle_id'])==pid and r['model']==m) for m in ['P65','P9A']]
        body.append([str(pid)]+[LABELS[r['outcome']] for r in pair])
        colors.append(['#f1f4f7']+[COLORS[r['outcome']]+'28' for r in pair])
    tab=ax.table(cellText=body,cellColours=colors,colLabels=['Particle ID','P6.5 only','P9-A'],
                 colWidths=[.18,.41,.41],cellLoc='center',bbox=[0,.06,1,.84])
    tab.auto_set_font_size(False);tab.set_fontsize(12)
    for (i,j),cell in tab.get_celld().items():
        cell.set_edgecolor('white')
        if i==0:cell.set_facecolor('#29475c');cell.set_text_props(color='white',weight='bold')
    ax.set_title('Same 12 particles | Same 2.0 mm/s FEM field',fontsize=17,pad=20)
    ax.text(.5,.96,'Completed: P6.5 = 5/12   |   P9-A = 5/12',ha='center',fontsize=13,transform=ax.transAxes)
    fig.text(.5,.03,'Identical births, radii, admission, time step, refinement and stop rules.',ha='center',fontsize=10)
    save(fig,'01_same12_p65_vs_p9a',['data/ab_outcomes.csv'])

def inlet2():
    rows=list(csv.DictReader((DATA/'inlet_escape_states.csv').open()))
    fig,axes=plt.subplots(2,3,figsize=(14,7.5),layout='constrained')
    plotrows=[]
    for i,summary in enumerate(S['inlet']):
        pid=summary['particle_id'];ax,bx,cx=axes[i];end=.001
        for model,color in [('P65','#236a9f'),('P9A',COLORS['INLET_ESCAPE'])]:
            trs=[r for r in rows if int(r['particle_id'])==pid and r['model']==model and r['record_kind']=='TRIAL_VELOCITY_EVALUATION' and r['trial_accepted']=='True' and float(r['time_s'])<=end]
            times=[];speeds=[]
            for r in trs:
                t=float(r['time_s']);t1=min(t+float(r['accepted_dt_s']),end)
                times.extend([t*1000,t1*1000]);speeds.extend([float(r['velocity_dot_inlet_inward_normal'])*1000]*2)
                plotrows.append(dict(particle_id=pid,panel='inward_velocity',model=model,t0_s=t,t1_s=t1,value_m_s=float(r['velocity_dot_inlet_inward_normal'])))
            ax.plot(times,speeds,color=color,lw=2,label='P6.5' if model=='P65' else 'P9-A')
            st=[r for r in rows if int(r['particle_id'])==pid and r['model']==model and r['record_kind']=='ACCEPTED_STATE' and float(r['time_s'])<=end]
            cx.plot([float(r['time_s'])*1000 for r in st],[float(r['gap_ratio']) for r in st],'.-',color=color,label='P6.5' if model=='P65' else 'P9-A',markersize=4)
        ax.hlines(summary['FEM_inward_velocity_m_s']*1000,0,.25,color='#555',ls='--',lw=2,label='FEM at P9-A start')
        ax.axhline(0,color='#999',lw=.7);ax.set_xlim(0,1)
        ax.text(.05,.18,f"P9-A: {summary['P9A_inward_velocity_m_s']*1000:.2e} mm/s\n(outward; exit after first step)",transform=ax.transAxes,color=COLORS['INLET_ESCAPE'],fontsize=10)
        ax.set_title(f'ID {pid} | Inward velocity');ax.set_ylabel('Velocity (mm/s)')
        values=[summary['bulk_t_speed'],summary['linear_shear_speed'],summary['target_tangential_speed']]
        for v,label,color,style in zip(values,[r'$|u_t|$',r'$|H s|$','P9-A target'],['#236a9f','#d3871b','#c33c49'],['-','--',':']):
            bx.plot([0,.25],[v*1000]*2,label=label,color=color,ls=style,lw=2.5)
            plotrows.append(dict(particle_id=pid,panel='closure_speeds',quantity=label,t0_s=0,t1_s=.00025,value_m_s=v))
        bx.set_yscale('log');bx.set_xlim(0,.25);bx.set_title('First accepted interval');bx.set_ylabel('Tangential speed (mm/s; log)')
        cx.set_title('Wall proximity');cx.set_ylabel('Wall gap / radius');cx.set_xlim(0,1)
        cx.text(.04,.92,'P9-A wall weight = 1',transform=cx.transAxes,va='top',color='#7a2630')
        for aa in [ax,bx,cx]:aa.set_xlabel('Elapsed time from birth (ms)');aa.grid(alpha=.18)
        if i==0:
            ax.legend(loc='upper right',fontsize=8);bx.legend(loc='center right',fontsize=9);cx.legend(fontsize=9)
    fig.suptitle('Inlet escapes | The velocity changes in the very first solve',fontsize=16)
    csvwrite(DATA/'figure02_plotted_intervals.csv',plotrows)
    save(fig,'02_inlet_escape_diagnosis',['data/inlet_escape_states.csv','data/inlet_summary.csv','data/figure02_plotted_intervals.csv'])

def geometry3():
    from particle_3d.particle81_simulation import environment
    env=environment();pid=3;m,a,t,states=load('P9A',pid);_,pa,_,_=load('P65',pid)
    first=t[0];center=a[0,1:4];radius=m['radius_m'];surf=env.boundaries['INLET']
    cap=np.asarray(surf.points)[surf.faces.reshape(-1,4)[:,1:]]
    wall=env.wall.triangles;scale=4*radius
    wi=np.flatnonzero(np.linalg.norm(wall.mean(axis=1)-center,axis=1)<scale)
    ci=np.flatnonzero(np.linalg.norm(cap.mean(axis=1)-center,axis=1)<scale)
    nearest=first['nearest_wall_triangle_id'];wi=np.unique(np.r_[wi,nearest])
    pshort=pa[pa[:,0]<=.001]
    geometry=dict(particle_id=pid,coordinate_role='ORIGINAL_FROZEN_MESH_SI_TRANSLATED_BY_BIRTH_CENTER_FOR_PLOT',
        origin_m=center.tolist(),wall_triangle_ids=wi.tolist(),wall_triangles_m=wall[wi].tolist(),
        inlet_cap_triangle_ids=ci.tolist(),inlet_cap_triangles_m=cap[ci].tolist(),nearest_wall_triangle_id=nearest,
        nearest_wall_triangle_m=wall[nearest].tolist(),radius_m=radius,P9A_positions_m=a[:,1:4].tolist(),
        P65_positions_m=pshort[:,1:4].tolist(),normal=first['wall_normal_xyz'],FEM_velocity_m_s=first['FEM_velocity_xyz'],
        shear_vector_s_inv=first['local_shear_vector_xyz'],arrow_role='NORMALIZED_DIRECTIONS; EACH_ARROW_LENGTH_1P4_RADII',
        cap_signed_distances_m=[S['inlet'][0]['distance_to_inlet_plane_m'],S['inlet'][0]['P9A_last_signed_distance_m']],
        P9A_times_s=a[:,0].tolist())
    dump(DATA/'inlet_geometry_closeup.json',geometry)
    fig=plt.figure(figsize=(13,6.4));ax=fig.add_subplot(121,projection='3d',computed_zorder=False);bx=fig.add_subplot(122)
    trans=lambda x:(np.asarray(x)-center)*1e6
    ax.add_collection3d(Poly3DCollection(trans(wall[wi]),facecolor='#7795aa',edgecolor='#627d91',linewidth=.35,alpha=.15,zorder=1))
    ax.add_collection3d(Poly3DCollection(trans(cap[ci]),facecolor='#d9e0e5',edgecolor='#697b88',linewidth=.5,alpha=.25,zorder=2))
    ax.add_collection3d(Poly3DCollection(trans(wall[[nearest]]),facecolor='#e4a33d',edgecolor='#9a6310',linewidth=1,alpha=.8,zorder=3))
    uu,vv=np.meshgrid(np.linspace(0,2*np.pi,28),np.linspace(0,np.pi,16))
    ax.plot_wireframe(radius*1e6*np.cos(uu)*np.sin(vv),radius*1e6*np.sin(uu)*np.sin(vv),radius*1e6*np.cos(vv),color='#b84c62',alpha=.3,lw=.45)
    p=trans(pshort[:,1:4]);ax.plot(*p.T,color='#236a9f',lw=2,marker='.',markersize=3,zorder=6)
    ax.scatter(0,0,0,color='#c33c49',s=40,depthshade=False,zorder=9)
    for key,col,label in [('wall_normal_xyz','#333','n'),('FEM_velocity_xyz','#248477','u'),('local_shear_vector_xyz','#a55dad','s')]:
        direction=np.array(first[key]);direction/=np.linalg.norm(direction)
        ax.quiver(0,0,0,*direction,length=1.4*radius*1e6,color=col,linewidth=2.5,arrow_length_ratio=.2,zorder=8)
        ax.text(*(direction*1.65*radius*1e6),label,color=col,fontsize=12,weight='bold',zorder=10)
    lim=scale*1e6
    ax.set(xlim=(-lim,lim),ylim=(-lim,lim),zlim=(-lim,lim),xlabel='x - birth x (um)',ylabel='y - birth y (um)',zlabel='z - birth z (um)')
    ax.set_box_aspect((1,1,1));ax.view_init(elev=24,azim=-53);ax.set_title('Actual inlet mesh and first positions',pad=20)
    bx.plot(a[:,0]*1000,np.array(geometry['cap_signed_distances_m'])*1e12,'o-',color='#c33c49',lw=2)
    bx.axhline(0,color='#444',ls='--',lw=1);bx.set_xlabel('Elapsed time from birth (ms)');bx.set_ylabel('Signed distance to actual cap triangle (pm)')
    bx.set_title('P9-A crosses the inlet cap');bx.grid(alpha=.2)
    bx.text(.05,.1,'Negative = outside\nP9-A displacement is too small\nto resolve at mesh scale.',transform=bx.transAxes,fontsize=11)
    handles=[Patch(color='#7795aa',alpha=.4,label='WALL triangles'),Patch(color='#c4ced6',label='Inlet cap'),Patch(color='#e4a33d',label='Nearest WALL triangle'),
        Line2D([0],[0],color='#236a9f',label='P6.5 first 1 ms'),Line2D([0],[0],marker='o',color='#c33c49',label='Birth / P9-A positions'),
        Line2D([0],[0],color='#333',label='Wall normal n'),Line2D([0],[0],color='#248477',label='FEM velocity u'),Line2D([0],[0],color='#a55dad',label='Shear s')]
    fig.legend(handles=handles,loc='lower center',ncol=4,fontsize=9,bbox_to_anchor=(.5,.015),frameon=False)
    fig.text(.5,.12,'Geometry to scale; arrows show normalized directions (not relative magnitudes).',ha='center',fontsize=9)
    fig.suptitle('ID 3 | Real inlet-rim geometry',fontsize=16);fig.subplots_adjust(bottom=.24,wspace=.34,top=.86)
    save(fig,'03_inlet_geometry_closeup',['data/inlet_geometry_closeup.json','data/inlet_summary.csv'])

def timeline4():
    fig,axes=plt.subplots(4,2,figsize=(12,11),layout='constrained');plotted=[]
    for row,summary in zip(axes,S['stall']):
        pid=summary['particle_id'];m,a,trials,_=load('P9A',pid);trials=trials[-100:];stop=a[-1,0]
        for accepted,marker,color in [(True,'o','#236a9f'),(False,'x','#c33c49')]:
            rs=[r for r in trials if r['trial_accepted']==accepted]
            row[0].scatter([(r['time_s']-stop)*1e6 for r in rs],[r['g_nf_m']*1e12 for r in rs],s=20,marker=marker,color=color,label='Accepted trial' if accepted else 'Rejected / stop trial')
        acc=[r for r in trials if r['trial_accepted']]
        row[1].plot([(r['time_s']-stop)*1e6 for r in acc],[r['accepted_dt_s']/r['requested_dt_s'] for r in acc],'o-',markersize=3,lw=1,color='#236a9f')
        row[1].set_yscale('log');row[0].axhline(0,color='#777',lw=.7)
        row[0].set_title(f'ID {pid} | Stop at {stop*1000:.6f} ms',loc='left')
        row[0].set_ylabel('Gap - handoff gap (pm)');row[1].set_ylabel('Accepted dt / nominal dt')
        row[1].text(.04,.91,f"{summary['P9A_last100']['triangle_switch_count']} triangle switches",transform=row[1].transAxes,va='top')
        for r in trials:
            if r['wall_triangle_changed']:
                for ax in row:ax.axvline((r['time_s']-stop)*1e6,color='#555',ls=':')
            plotted.append(dict(particle_id=pid,trial_index=r['trial_index'],time_s=r['time_s'],
                time_before_stop_us=(r['time_s']-stop)*1e6,g_nf_pm=r['g_nf_m']*1e12,
                accepted=r['trial_accepted'],accepted_dt_over_nominal_dt=r['accepted_dt_s']/r['requested_dt_s'],wall_triangle_changed=r['wall_triangle_changed']))
        for ax in row:ax.set_xlabel('Time relative to stop (us)');ax.grid(alpha=.18);ax.ticklabel_format(axis='x',style='plain',useOffset=False)
    axes[0,0].legend(fontsize=8,loc='upper right')
    fig.suptitle('Handoff stalls | Last 100 P9-A trials for each particle',fontsize=16)
    csvwrite(DATA/'figure04_plotted_trials.csv',plotted)
    save(fig,'04_handoff_stall_timeline',['data/handoff_stall_trials.csv','data/figure04_plotted_trials.csv','data/analysis_summary.json'])

def summary5():
    fig,ax=plt.subplots(figsize=(15,5.2));ax.axis('off');body=[]
    for r in S['stall']:
        z=r['P9A_last100']
        body.append([str(r['particle_id']),LABELS[r['P65_outcome']],LABELS[r['P9A_outcome']],
            LABELS[r['checkpoint_P65_outcome']],str(z['triangle_switch_count']),f"{z['maximum_normal_change_deg']:.4f}",
            f"{r['minimum_accepted_g_nf_m']:.2e}",f"{z['condition_min']:.1f} - {z['condition_max']:.1f}"])
    tab=ax.table(cellText=body,colLabels=['ID','P6.5\nfrom birth','P9-A\nfrom birth','P6.5 from same\nP9-A checkpoint','Triangle\nswitches','Max normal\nchange (deg)','Min accepted\ng_nf (m)','Scaled matrix\ncondition'],
        cellLoc='center',colWidths=[.045,.155,.155,.185,.08,.11,.13,.14],bbox=[0,.24,1,.52])
    tab.auto_set_font_size(False);tab.set_fontsize(10)
    for (i,j),cell in tab.get_celld().items():
        cell.set_edgecolor('white')
        if i==0:cell.set_facecolor('#29475c');cell.set_text_props(color='white',weight='bold')
        else:cell.set_facecolor('#edf2f6' if i%2 else '#f7f9fb')
    ax.set_title('Four stalls | Endpoint handoff rejection dominates',fontsize=17,pad=15)
    ax.text(.5,.86,'Dominant rejection for all four: HANDOFF_ENDPOINT_BELOW_LOWER',ha='center',transform=ax.transAxes,fontsize=12,color='#974b24')
    ax.text(0,.12,'Switches, normal changes and resistance condition ranges: last 100 P9-A trials. Minimum g_nf: all accepted states.\nSmall negative g_nf values lie within the original geometric roundoff budget (~2.01e-17 m).',transform=ax.transAxes,fontsize=10)
    ax.text(0,.015,'ID 7 also has a separate contact-system rank issue: condition = 1.19e13 for two constraints on the same edge.',transform=ax.transAxes,fontsize=10,color='#974b24')
    save(fig,'05_handoff_failure_reason_summary',['data/analysis_summary.json','data/handoff_stall_trials.csv','data/ab_outcomes.csv'])

def consistency6():
    rows=list(csv.DictReader((DATA/'closure_consistency_states.csv').open()));fig,ax=plt.subplots(figsize=(9,7))
    for outcome in LABELS:
        rs=[r for r in rows if r['outcome']==outcome]
        ax.scatter([float(r['bulk_t_speed'])*1000 for r in rs],[float(r['linear_shear_speed'])*1000 for r in rs],
            color=COLORS[outcome],s=26 if outcome in ['COMPLETED','SUBDIVISION_STOP'] else 55,
            alpha=.7,edgecolor='white',linewidth=.3,label=f'{LABELS[outcome]} ({len(rs)} states)',zorder=3)
    for r in S['inlet']:
        ax.annotate(f"ID {r['particle_id']}",(r['bulk_t_speed']*1000,r['linear_shear_speed']*1000),
            xytext=(-55,-17 if r['particle_id']==3 else 17),textcoords='offset points',fontsize=10,
            arrowprops=dict(arrowstyle='-',color='#777',lw=.6))
    xmax=max(float(r['bulk_t_speed'])*1000 for r in rows)*1.2
    ymin=min(float(r['linear_shear_speed'])*1000 for r in rows)/2
    ymax=max(xmax,max(float(r['linear_shear_speed'])*1000 for r in rows)*1.2)
    ax.plot([.5,xmax],[.5,xmax],'--',color='#444',lw=1,label='y = x')
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(.5,xmax);ax.set_ylim(ymin,ymax)
    ax.set_xticks([.5,1,2,4,6],labels=['0.5','1','2','4','6'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    assert all(.5<float(r['bulk_t_speed'])*1000<xmax and ymin<float(r['linear_shear_speed'])*1000<ymax for r in rows)
    dump(DATA/'figure06_axis_audit.json',dict(plotted_count=len(rows),all_points_within_axes=True,
        x_limits_mm_s=[.5,xmax],y_limits_mm_s=[ymin,ymax],log_axes=True))
    ax.set_xlabel(r'FEM tangential speed $|u_t|$ (mm/s)');ax.set_ylabel(r'Local linear-shear prediction $|H s|$ (mm/s)')
    ax.set_title('Bulk velocity vs local planar shear | P9-A states',pad=15,fontsize=15)
    ax.legend(loc='center left',fontsize=9,frameon=True,facecolor='white',edgecolor='#ddd')
    ax.grid(which='major',alpha=.18)
    fig.text(.5,.015,'Up to 25 accepted evaluation states per particle; colors denote the full-trajectory outcome. No fitted threshold.',ha='center',fontsize=9)
    fig.subplots_adjust(bottom=.12)
    save(fig,'06_bulk_vs_shear_consistency',['data/closure_consistency_states.csv','data/analysis_summary.json'])

if __name__=='__main__':
    table1();inlet2();geometry3();timeline4();summary5();consistency6()
    dump(DATA/'figure_sources.json',MANIFEST)
