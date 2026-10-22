"""Five reproducible 300 dpi figures from unsmoothed numerical samples."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pyvista as pv
from analyze import quantities

LABELS=['M0 original','Fixed PLC medium','Fixed PLC fine']
COLORS=['#252525','#1273ab','#cf631e']
EVENT_COLORS=['#b2182b','#2166ac','#1b7837','#762a83']
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,'legend.fontsize':9,'figure.facecolor':'white','savefig.facecolor':'white'})
P=read(HERE/'EVENT_FIELD_PROTOCOL.json');CONFIG=read(HERE/'SPATIAL_CONFIGURATION.json');MET=read(HERE/'FIELD_METRICS.json');FIG={}

def save(fig,name,description):
    fig.savefig(HERE/name,dpi=300,bbox_inches='tight')
    plt.close(fig);FIG[name]=description;print('FIGURE',name,flush=True)

def frame(ax,xyz):
    lo=xyz.min(0);hi=xyz.max(0);c=(lo+hi)/2;r=(hi-lo).max()/2*1.06
    ax.set(xlim=(c[0]-r,c[0]+r),ylim=(c[1]-r,c[1]+r),zlim=(c[2]-r,c[2]+r))
    ax.set_box_aspect((1,1,1));ax.view_init(elev=22,azim=-65)
    ax.set_xlabel('x [mm]');ax.set_ylabel('y [mm]');ax.set_zlabel('z [mm]')

def surface(ax,grid,origin,alpha=.1):
    if not isinstance(grid,pv.PolyData):grid=grid.extract_surface(algorithm='dataset_surface')
    tri=grid.points[grid.faces.reshape(-1,4)[:,1:]]
    c=Poly3DCollection((tri-origin)*1000,facecolor='#9ba6af',edgecolor='none',alpha=alpha,rasterized=True)
    ax.add_collection3d(c)

def event_regions():
    fig=plt.figure(figsize=(17,10));gs=fig.add_gridspec(2,3,width_ratios=[1.3,1,1])
    ax=fig.add_subplot(gs[:,0],projection='3d');wall=pv.read(HERE/'vtk/original_wall.vtp');surface(ax,wall,np.zeros(3),.14)
    for pid,c in zip(PIDS,EVENT_COLORS):
        z=np.load(OLD/f'original_{pid}.npz');xyz=z['x_m']*1000
        ax.plot(*xyz.T,color=c,lw=2,label=f'ID {pid}')
        point=np.array(P['events'][str(pid)]['A']['x_m'])*1000;ax.scatter(*point,color=c,s=25);ax.text(*point,str(pid),color=c)
    frame(ax,wall.points*1000);ax.set_title('Original BraVa geometry\nHistorical accepted paths; no new integration');ax.legend(loc='upper left')
    for j,(pid,c) in enumerate(zip(PIDS,EVENT_COLORS)):
        ax=fig.add_subplot(gs[j//2,1+j%2],projection='3d');e=P['events'][str(pid)];origin=np.array(e['A']['x_m'])
        roi=pv.read(HERE/f'vtk/original_roi_{pid}.vtu');surface(ax,roi,origin,.13)
        z=np.load(OLD/f'original_{pid}.npz');x=(z['x_m']-origin)*1000
        ax.plot(*x.T,c=c,lw=1.5)
        for label,marker in [('A','o'),('B','s')]:
            point=(np.array(e[label]['x_m'])-origin)*1000
            ax.scatter(*point,c=c,s=25,marker=marker);ax.text(*point,label)
            u=np.array(MET[str(pid)]['representatives']['M0'][label]['velocity_m_s'])*1000*.005
            ax.quiver(*point,*u,color='black',arrow_length_ratio=.25)
        plane=CONFIG[str(pid)]['plane'];o=np.array(plane['origin_m']);e1=np.array(plane['e1']);e2=np.array(plane['e2']);lo,hi=np.array(plane['bounds_m'])
        corners=np.array([o+s*e1+t*e2 for s,t in [(lo[0],lo[1]),(hi[0],lo[1]),(hi[0],hi[1]),(lo[0],hi[1])]])
        ax.add_collection3d(Poly3DCollection([(corners-origin)*1000],facecolor='gold',edgecolor='#a87800',alpha=.35))
        frame(ax,(roi.points-origin)*1000);ax.set_title(f'ID {pid}: path + two neighbor layers\nLocal coordinates relative to A; arrows = u x 5 ms')
    fig.suptitle('Event locations in the original steady field — fixed world coordinates',fontsize=16)
    save(fig,'01_EVENT_REGIONS.png',dict(camera={'elev':22,'azim':-65,'equal_3D_aspect':True},source='original_wall.vtp, original_roi_*.vtu, original_*.npz; original representative u',local_origins='each original A',arrows_time_s=.005,plane='SPATIAL_CONFIGURATION.json'))

def plane_overlay(ax,pid):
    e=P['events'][str(pid)];pl=CONFIG[str(pid)]['plane'];o=np.array(pl['origin_m']);e1=np.array(pl['e1']);e2=np.array(pl['e2']);basis=np.stack([e1,e2],1)
    roi=pv.read(HERE/f'vtk/original_roi_{pid}.vtu');cut=roi.slice(normal=pl['normal'],origin=o)
    edge=cut.extract_all_edges();lines=edge.lines;offset=0
    while offset<len(lines):
        n=lines[offset];ids=lines[offset+1:offset+n+1];p=(edge.points[ids]-o)@basis*1e6
        ax.plot(p[:,0],p[:,1],c='.3',lw=.35,alpha=.6);offset+=n+1
    pt=(np.array(e['A']['x_m'])-o)@basis*1e6
    ax.add_patch(Circle(pt,e['radius_m']*1e6,facecolor='none',edgecolor='#ef28e3',lw=1.4))
    ax.scatter(*pt,c='#ef28e3',s=10,zorder=10)
    ax.annotate('A',pt,xytext=(5,5),textcoords='offset points',color='#b600a9')
    if pid==1189:
        wall=pv.read(HERE/'vtk/original_wall.vtp');tri=wall.points[wall.faces.reshape(-1,4)[33200,1:]]
        p=(np.r_[tri,tri[:1]]-o)@basis*1e6
        ax.plot(*p.T,c='limegreen',lw=1.5)
    lo,hi=np.array(pl['bounds_m'])*1e6
    ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),aspect='equal')
    ax.set_xlabel('tangent coordinate [um]');ax.set_ylabel('inward normal coordinate [um]')

def near_wall():
    fig=plt.figure(figsize=(18,15));gs=fig.add_gridspec(4,3,height_ratios=[1,1,.8,.8],hspace=.5,wspace=.3)
    pid=1189;plane=np.load(HERE/f'plane_{pid}.npz');e=P['events'][str(pid)];n=np.array(e['plane']['e2']);tangent=np.array(e['plane']['e1'])
    fields=[np.load(HERE/f'{name}_plane_{pid}.npz')['flow'] for name in NAMES]
    vmax=max(np.nanmax(norm(f[:,:3])) for f in fields)*1000;unmax=max(np.nanmax(abs(f[:,:3]@n)) for f in fields)*1000
    x=plane['axis1_m']*1e6;y=plane['axis2_m']*1e6;axes=[[],[]]
    for j,(name,f) in enumerate(zip(NAMES,fields)):
        for row in [0,1]:
            ax=fig.add_subplot(gs[row,j]);axes[row].append(ax)
            val=norm(f[:,:3])*1000 if row==0 else f[:,:3]@n*1000
            im=ax.pcolormesh(x,y,val.reshape(x.shape),shading='nearest',cmap='viridis' if row==0 else 'RdBu_r',vmin=0 if row==0 else -unmax,vmax=vmax if row==0 else unmax,rasterized=True)
            uf=(f[:,:3]@tangent*1000).reshape(x.shape);vf=(f[:,:3]@n*1000).reshape(x.shape)
            ax.quiver(x[::16,::16],y[::16,::16],uf[::16,::16],vf[::16,::16],angles='xy',scale_units='xy',scale=2.,color='white' if row==0 else 'black',width=.003)
            plane_overlay(ax,pid);ax.set_title(f'1189: {LABELS[j]} | '+('speed' if row==0 else 'signed inward component'))
        fig.colorbar(im,ax=axes[1][-1],shrink=.68,label='u dot n [mm/s]')
    # One shared speed bar; the normal panels use the same symmetric limits.
    import matplotlib.cm as cm
    fig.colorbar(cm.ScalarMappable(norm=plt.Normalize(0,vmax),cmap='viridis'),ax=axes[0],shrink=.72,label='speed [mm/s]')
    for col,pid in enumerate([1189,1329,642]):
        coord=np.load(HERE/f'profiles_{pid}.npz');wall=np.load(HERE/f'wall_profiles_{pid}.npz');e=P['events'][str(pid)]
        for row,k in [(2,0),(3,1)]:
            ax=fig.add_subplot(gs[row,col]);mask=(coord['profile_index']==k)&wall['original_normal_ray_valid'];d=coord['d_m']*1e6
            n=np.array(e['normal_profiles']['normals'][k])
            for name,color,label in zip(NAMES,COLORS,LABELS):
                f=np.load(HERE/f'{name}_profiles_{pid}.npz')['flow'];u=f[:,:3].copy();u[coord['d_m']==0]=0
                un=u@n;ut=norm(u-un[:,None]*n)
                ax.plot(d[mask],un[mask]*1000,c=color,label=label+' normal')
                ax.plot(d[mask],ut[mask]*1000,c=color,ls='--',alpha=.9,label=label+' tangent')
            ax.axhline(0,c='.6',lw=.5);ax.axvline(e['radius_m']*1e6,c='#b600a9',lw=.7,label='sphere radius')
            ax.set_title(f'ID {pid}, normal ray from '+('A (slow)' if k==0 else 'B (fast)'))
            ax.set_xlabel('d: original wall distance [um]');ax.set_ylabel('components [mm/s]');ax.grid(alpha=.2)
    handles,labels=ax.get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.015),ncol=4,fontsize=8)
    fig.suptitle('Original near-wall P1 input on matched geometry\nGreen: original facet 33200; wireframe: original cells including 388769; magenta: actual sphere center/radius\nArrows: in-plane velocity, 1 mm/s = 0.5 um of arrow length; profile gap h = d - radius',fontsize=13,y=.99)
    save(fig,'02_NEAR_WALL_FLOW.png',dict(speed_limits_mm_s=[0,vmax],normal_limits_mm_s=[-unmax,unmax],saturated_points=0,plane='SPATIAL_CONFIGURATION.json',arrows_scale_mm_s_per_um=2.,profile_trace_at_wall='Verified no-slip u=0 only; G not assigned by trace',profile_stop='original_normal_ray_valid mask',source='*_plane_1189.npz, *_profiles_*.npz'))

def fixed_paths():
    fig,axes=plt.subplots(4,4,figsize=(19,13),constrained_layout=True)
    for row,pid in enumerate(PIDS):
        path=np.load(HERE/f'path_{pid}.npz');ix=path['quadrature_order']==1;e=P['events'][str(pid)]
        n=np.load(HERE/f'wall_path_{pid}.npz')['normal']
        if pid==822:n=np.broadcast_to(e['plane']['e1'],n.shape)
        qs={name:quantities(np.load(HERE/f'{name}_path_{pid}.npz')['flow'],n) for name in NAMES}
        for phase in [0,1,2]:
            ax=axes[row,phase];mask=ix&(path['phase']==phase);s=path['arc_m'][mask]*1e6
            for name,color,label in zip(NAMES,COLORS,LABELS):
                ax.plot(s,qs[name]['speed'][mask]*1000,c=color,lw=1,label=label)
                ax.plot(s,qs[name]['un'][mask]*1000,c=color,ls='--',lw=.85)
            ax.set_title(f'{pid} | {e["phase_names"][phase]}');ax.set_xlabel('path arc length [um]');ax.set_ylabel('speed / component [mm/s]');ax.grid(alpha=.2)
            if row==0 and phase==0:ax.legend()
        ax=axes[row,3];s=path['arc_m'][ix]*1e6
        for name,color in zip(NAMES[1:],COLORS[1:]):
            ax.plot(s,norm(qs[name]['u'][ix]-qs['M0']['u'][ix])*1000,c=color,lw=1,label=name.replace('M1_fixedplc_','')+' vs M0')
        ax.set_title(f'{pid} | vector difference from M0');ax.set_xlabel('path arc length [um]');ax.set_ylabel('|delta u| [mm/s]');ax.grid(alpha=.2)
        if row==0:ax.legend()
    fig.suptitle('Same historical spatial path in three independently solved fields — linear axes\nSolid: |u|; dashed: inward component (1189/1329/642), fixed path-direction component (822). Phase scales differ; no new particle velocities.',fontsize=14)
    save(fig,'03_FIXED_PATH_FIELD_COMPARISON.png',dict(source='Common union-span midpoints; *_path_*.npz',phases='Frozen original phase boundaries',units='um, mm/s',axes='Independent linear vertical scales for each phase; all fields share a panel scale'))

def rise822():
    pid=822;fig=plt.figure(figsize=(17,10));gs=fig.add_gridspec(2,3)
    e=P['events'][str(pid)];a=np.array(e['A']['x_m']);b=np.array(e['B']['x_m']);r=MET[str(pid)]['spatial_rise'];parts=read(HERE/'PARTITIONS.json')
    ax=fig.add_subplot(gs[:,0],projection='3d');roi=pv.read(HERE/'vtk/original_roi_822.vtu')
    ids=[int(z[2]) for z in r['M0']['own_cells']];mask=np.isin(roi.cell_data['original_canonical_cell'],ids);cells=roi.extract_cells(mask)
    surface(ax,cells,a,.2);edge=cells.extract_all_edges();ls=edge.lines;off=0
    while off<len(ls):
        k=ls[off];xyz=(edge.points[ls[off+1:off+k+1]]-a)*1000;ax.plot(*xyz.T,c='.3',lw=.7);off+=k+1
    ax.plot(*np.array([np.zeros(3),(b-a)*1000]).T,c='#b2182b',lw=3);ax.text(0,0,0,'A');ax.text(*((b-a)*1000),'B')
    frame(ax,(cells.points-a)*1000);ax.set_title('Original four tetrahedra\nFixed A-to-B segment, 168.6564 um')
    coord=np.load(HERE/'rise_822.npz');s=coord['arc_m']*1e6
    for row in [0]:
        ax=fig.add_subplot(gs[row,1:])
        for k,(name,color,label) in enumerate(zip(NAMES,COLORS,LABELS)):
            f=np.load(HERE/f'{name}_rise_822.npz')['flow'];val=norm(f[:,:3])*1000 if row==0 else f[:,:3]@np.array(e['plane']['e1'])*1000
            ax.plot(s,val,c=color,label=label,lw=1.8)
            for span in r[name]['own_cells'][:-1]:ax.axvline(span[1]*s[-1],c=color,ls=':',alpha=.4,lw=.8)
        if row==0:
            for j,level in enumerate(r['M0']['levels']):ax.axhline(level['level_m_s']*1000,c='.5',ls='--',label=f'original {10 if j==0 else 90}% level')
            text='\n'.join(f'{label}: gain {MET[str(pid)]["representatives"][name]["scalar_speed_increase_mm_s"]:.3f} mm/s; L10-90 '+(f'{r[name]["length_10_90_m"]*1e6:.2f} um' if r[name]['length_10_90_m'] is not None else 'not reached / ambiguous') for name,label in zip(NAMES,LABELS))
            ax.text(.02,.97,text,transform=ax.transAxes,va='top',fontsize=10,bbox={'facecolor':'white','alpha':.85,'edgecolor':'.8'})
        ax.set(xlabel='distance from original A [um]',ylabel='speed [mm/s]' if row==0 else 'fixed path component [mm/s]');ax.grid(alpha=.2);ax.legend(loc='lower right')
    sub=gs[1,1:].subgridspec(1,3,wspace=.32)
    plane=np.load(HERE/'plane_822.npz');xx=plane['axis1_m']*1e6;yy=plane['axis2_m']*1e6
    ff=[np.load(HERE/f'{name}_plane_822.npz')['flow'] for name in NAMES]
    vmin=min(np.nanmin(norm(f[:,:3])) for f in ff)*1000;vmax=max(np.nanmax(norm(f[:,:3])) for f in ff)*1000
    pa=[]
    for j,(name,label,f) in enumerate(zip(NAMES,LABELS,ff)):
        ax=fig.add_subplot(sub[0,j]);pa.append(ax)
        im=ax.pcolormesh(xx,yy,(norm(f[:,:3])*1000).reshape(xx.shape),shading='nearest',vmin=vmin,vmax=vmax,cmap='viridis',rasterized=True)
        ax.scatter(0,0,c='#ef28e3',s=15);ax.annotate('A',(0,0),xytext=(5,5),textcoords='offset points',color='#ef28e3')
        ax.set_title(label+'\nSame independent plane');ax.set_aspect('equal');ax.set_xlabel('path-frame axis [um]')
        if j==0:ax.set_ylabel('transverse axis [um]')
    fig.colorbar(im,ax=pa,shrink=.65,label='speed [mm/s]')
    fig.suptitle('ID 822: spatial field validation, not a new timestep experiment\nDotted lines: each field’s actual cell interfaces; identical original 10–90% levels for all fields',fontsize=15)
    save(fig,'04_ID822_SPATIAL_RISE.png',dict(source='rise_822.npz, *_rise_822.npz, *_plane_822.npz, FIELD_METRICS.json',camera={'elev':22,'azim':-65,'equal_3D_aspect':True},origin_m=a.tolist(),fixed_levels_m_s=[x['level_m_s'] for x in r['M0']['levels']],no_particle_dt_change=True,independent_plane_speed_limits_mm_s=[vmin,vmax],saturated_points=0))

def balances():
    fig,axes=plt.subplots(3,4,figsize=(19,12),constrained_layout=True)
    balance={name:read(HERE/f'{name}_BALANCE.json') for name in NAMES}
    for col,pid in enumerate(PIDS):
        for name,color,label in zip(NAMES,COLORS,LABELS):
            b=[r for r in balance[name] if r['event']==pid];x=[r['volume_m3']*1e18 for r in b]
            q=[r['levels'][-1]['net_flux_m3_s'] for r in b];through=[r['levels'][-1]['throughput_m3_s'] for r in b]
            err=[r['final_face_L1_refinement_difference_m3_s'] for r in b]
            axes[0,col].errorbar(x,np.array(q)*1e12,yerr=np.array(err)*1e12,c=color,marker='o',label=label)
            axes[1,col].plot(x,[qq/ss if ss else np.nan for qq,ss in zip(q,through)],c=color,marker='o')
        axes[0,col].set_title(f'{pid}: same closed surfaces');axes[0,col].set_ylabel('Qnet [1e-12 m3/s]');axes[0,col].set_xlabel('CV volume [um3]')
        axes[1,col].set_ylabel('Qnet / (0.5 sum |Qface|)');axes[1,col].set_xlabel('CV volume [um3]')
        for row in [0,1]:axes[row,col].axhline(0,c='.5',lw=.6);axes[row,col].grid(alpha=.2);axes[row,col].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
        ax=axes[2,col];ax.axis('off');e=MET[str(pid)];text=f'ID {pid} | medium vs fine\nArc-weighted RMS differences\n'
        for label in e['phase_names']:
            result=e['phase_comparisons'][NAMES[1]+'__'+NAMES[2]][label]['arc_weight_m'];m=result['metrics']
            text+=f'\n{label}:\n  u {m["u"]["RMS_difference"]*1000:.4g} mm/s\n  direction {m["direction_deg"]["RMS_difference"]:.3g} deg\n  G {m["G"]["RMS_difference"]:.4g} s^-1\n  coverage {result["valid_points"]}/{result["total_points"]}\n'
        text+='\nReference boundary / BC / steady: admitted\nNo physical validation; no convergence PASS'
        ax.text(0,1,text,va='top',fontsize=10)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Local continuity residual and resolution evidence\nALL CV boundary faces included. Error bars: final face-L1 quadrature change (not a physical uncertainty bound).\nExact div-volume identity checks arithmetic only; nonzero Qnet is the local continuity residual.',fontsize=14)
    save(fig,'05_LOCAL_BALANCE_AND_EVIDENCE.png',dict(source='*_BALANCE.json, FIELD_METRICS.json',normalization='0.5 sum abs per-original-face flux; no epsilon',flux_units='m3/s; panel multiplier 1e12',quadrature_levels=[0,1,2,3,4,5,6]))

def main():
    if len(sys.argv)>1:
        FIG.update(read(HERE/'FIGURE_METADATA.json'))
        functions=dict(event_regions=event_regions,near_wall=near_wall,fixed_paths=fixed_paths,rise822=rise822,balances=balances)
        for key in sys.argv[1:]:functions[key]()
    else:
        event_regions();near_wall();fixed_paths();rise822();balances()
    put('FIGURE_METADATA.json',FIG)
    links='\n'.join(f'<section><h2>{name}</h2><a href="{name}"><img src="{name}" style="width:100%;height:auto"></a></section>' for name in FIG)
    (HERE/'OPEN_RESULTS.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>BraVa 事件背景场核查</title><body style="max-width:1500px;margin:2em auto;font-family:sans-serif"><h1>BraVa 事件区域背景流场数值核查</h1><p>既有独立 CFD 共点比较；没有新 CFD 或微泡积分。</p><p><a href="FIELD_RELIABILITY_REPORT_ZH.md">中文报告</a> · <a href="FIELD_RELIABILITY_RESULTS.json">机器结果</a> · <a href="FIELD_COMPARABILITY.json">场准入</a> · <a href="REPRODUCE.md">复现与 ParaView 数据说明</a></p>'+links+'</body></html>')

if __name__=='__main__':main()
