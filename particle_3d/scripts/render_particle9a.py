"""Academic plots and an honest saved-state smoke replay (no new integration)."""
from pathlib import Path
import os,sys,json,csv,math
os.environ.setdefault('PYVISTA_OFF_SCREEN','true');os.environ.setdefault('LP_NUM_THREADS','2');os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'particle_3d/src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyvista as pv
import imageio_ffmpeg
from scipy.spatial.transform import Rotation,Slerp
from particle_3d.particle81_simulation import dump
from particle_3d.particle9a_provenance import FEM,sha256
REPORT=ROOT/'particle_3d/reports/particle9a_2mmps';OUT=ROOT/'particle_3d/outputs/particle9a_2mmps'
FIG=REPORT/'figures';ANI=REPORT/'animations'
for folder in [FIG,ANI]:folder.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.linewidth':.8,'savefig.facecolor':'white','figure.facecolor':'white','pdf.fonttype':42})
BLUE='#246da8';GOLD='#c47b10';RED='#b64248';GREEN='#218b83'


def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=230,bbox_inches='tight');fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');plt.close(fig)


def load_csv(name):
    with (REPORT/'data'/name).open() as f:return list(csv.DictReader(f))


def numeric(rows,key):return np.array([float(r[key]) for r in rows])


def scientific_plots():
    flow=json.loads((REPORT/'data/new_flow_flux.json').read_text())
    fig,ax=plt.subplots(1,2,figsize=(11,4.4),gridspec_kw={'width_ratios':[1,1.2]});fig.subplots_adjust(wspace=.35,bottom=.19,top=.86)
    ax[0].axis('off');ax[0].text(0,.99,'(a) Active background flow',weight='bold',transform=ax[0].transAxes)
    ax[0].text(0,.76,'2.0 mm/s',fontsize=31,color=BLUE,weight='bold',transform=ax[0].transAxes)
    ax[0].text(0,.58,'CURRENT INPUT  |  Solved steady FEM',fontsize=11,transform=ax[0].transAxes)
    ax[0].text(0,.39,f'Inlet Q = {flow["inlet_Q_uL_min"]:.6g} µL/min',transform=ax[0].transAxes)
    ax[0].text(0,.17,'Historical input: 0.352841 mm/s\nArchived; excluded from current production',color='#666666',fontsize=10,transform=ax[0].transAxes)
    roles=sorted(flow['outlet_fractions']);y=[100*flow['outlet_fractions'][r] for r in roles]
    ax[1].barh([r.replace('OUTLET_','Outlet ') for r in roles],y,color=[GREEN,BLUE,GOLD],height=.55);ax[1].invert_yaxis()
    for i,(r,v) in enumerate(zip(roles,y)):
        ax[1].text(v+1,i,f'{v:.2f}%',va='center',fontsize=11)
    ax[1].set(xlim=(0,105),xlabel='Flow relative to inlet (%)',title='(b) New outlet flow split');ax[1].grid(axis='x',alpha=.15);ax[1].set_axisbelow(True)
    fig.suptitle('New frozen flow | 2.0 mm/s',fontsize=16,weight='bold');save(fig,'01_current_flow')
    with (REPORT/'data/figure01_flow.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['boundary','Q_m3_s','Q_uL_min','fraction_of_inlet'])
        w.writerow(['INLET',flow['integrated_inlet_Q_m3_s'],flow['inlet_Q_uL_min'],1.])
        for r in roles:w.writerow([r,flow['outlet_flows_m3_s'][r],flow['outlet_flows_m3_s'][r]*60e9,flow['outlet_fractions'][r]])
    rows=load_csv('planar_comparison.csv');x=numeric(rows,'xi')
    fig,ax=plt.subplots(figsize=(8.6,4.8));fig.subplots_adjust(bottom=.22)
    for key,label,color in [('m_tt',r'$m_{tt}$ — tangential translation',BLUE),('m_rr',r'$m_{rr}$ — wall-parallel rotation',GREEN),('m_cross',r'$m_{cross}$ — reciprocal coupling',GOLD)]:ax.semilogx(x,numeric(rows,key),label=label,color=color,lw=2)
    ax.axvline(.1,color='#999999',ls=':',lw=1);ax.axvline(1,color='#999999',ls=':',lw=1)
    ax.set(xlabel=r'Gap ratio $h/a$',ylabel='Dimensionless effective mobility',title='Single planar wall | Sliding and rotation are coupled',ylim=(-.03,1.06));ax.legend(frameon=False,loc='upper left');ax.grid(alpha=.16)
    fig.text(.12,.035,'Near the wall, sliding and rotation change and influence each other.\nReciprocal projection uses the rounding range of the original coefficients.',fontsize=10,color='#555555');save(fig,'02_tangential_mobility')
    fig,ax=plt.subplots(1,3,figsize=(13.4,4.1));fig.subplots_adjust(wspace=.38,bottom=.25,top=.80)
    ax[0].semilogx(x,numeric(rows,'p65_Vt_mm_s'),color='#888888',ls='--',label='P6.5');ax[0].semilogx(x,numeric(rows,'p9a_Vt_mm_s'),color=BLUE,label='P9-A')
    ax[0].set(ylabel='Tangential velocity (mm/s)',title='(a) Sliding');ax[0].legend(frameon=False)
    ax[1].semilogx(x,numeric(rows,'p65_omega_y_s_inv'),color='#888888',ls='--',label='P6.5');ax[1].semilogx(x,numeric(rows,'p9a_omega_y_s_inv'),color=GOLD,label='P9-A')
    ax[1].set(ylabel=r'Angular velocity $\Omega_y$ (rad/s)',title='(b) Rotation');ax[1].legend(frameon=False)
    ax[2].semilogx(x,numeric(rows,'normal_difference_m_s'),color=GREEN);ax[2].set(ylabel='Normal velocity difference (m/s)',title='(c) P9-A minus P6.5',ylim=(-1e-12,1e-12));ax[2].text(.06,.82,'Maximum: 0 m/s',transform=ax[2].transAxes)
    for a in ax:a.set_xlabel(r'Gap ratio $h/a$');a.grid(alpha=.15)
    fig.suptitle('Planar shear case | Normal protection is unchanged',fontsize=15)
    fig.text(.07,.025,'Sphere radius: 1 µm; shear rate: 500 s⁻¹. Rotation sign follows the frozen 2D reference;\nnear-wall rotation changes sign relative to the bulk vorticity in this prescribed closure.',fontsize=10,color='#555555');save(fig,'03_p65_vs_p9a')
    summary=json.loads((OUT/'provenance/smoke_summary.json').read_text());r=summary['rows']
    travels=np.array([q['residence_time_s'] for q in r if q['completed']]);gaps=np.array([q['minimum_gap_m'] for q in r])*1e9
    fig,ax=plt.subplots(1,3,figsize=(13.4,4.3));fig.subplots_adjust(wspace=.38,top=.74,bottom=.27)
    ax[0].hist(travels*1000,bins=[20,30,40,50,60,70],color=BLUE,edgecolor='white');ax[0].set(xlabel='Travel time (ms)',ylabel='Completed tracks',title='(a) Completed only (n = 5)')
    labels=['Outlet 01','Outlet 02','Outlet 03','No exit'];counts=[0,5,0,7]
    ax[1].bar(np.arange(4),counts,color=[GREEN,BLUE,GOLD,RED]);ax[1].set_xticks(np.arange(4),labels,rotation=25,ha='right');ax[1].set(ylabel='Tracks',title='(b) All smoke outcomes (n = 12)',ylim=(0,8))
    for i,c in enumerate(counts):ax[1].text(i,c+.15,f'{c} ({100*c/12:.1f}%)',ha='center',fontsize=9)
    ax[2].hist(gaps,bins=np.linspace(0,100,11),color=GREEN,edgecolor='white');ax[2].set(xlabel='Minimum WALL gap (nm)',ylabel='Tracks',title='(c) All recorded tracks (n = 12)');ax[2].axvline(2,color=RED,ls='--',lw=1,label='Handoff: 2 nm');ax[2].legend(frameon=False,fontsize=9)
    fig.suptitle('2.0 mm/s + P9-A | Smoke batch failed the production gate',fontsize=16,color=RED)
    fig.text(.07,.075,'12 integrated: 5 outlet exits · 2 inlet escapes · 4 subdivision progress stops · 1 time-limit stop\n500-track production was not launched. No observed WALL penetration; upstream cap escape is a separate failure.',fontsize=11)
    save(fig,'05_dataset_status')
    with (REPORT/'data/figure05_smoke_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=r[0].keys());w.writeheader();w.writerows(r)


def replay_visuals():
    manifest=json.loads((FEM/'frozen_reference/boundary_manifest.json').read_text());wall=pv.read(FEM/manifest['boundaries']['WALL']['path'])
    origin=(wall.points.min(0)+wall.points.max(0))/2;wall.points=(wall.points-origin)*1e6
    records=[]
    for path in sorted((OUT/'trajectories').glob('mb_*.json')):
        if '.receipt.' in path.name:continue
        m=json.loads(path.read_text());s=np.load(OUT/m['samples_path'])['samples']
        if sha256(OUT/m['samples_path'])!=m['samples_sha256']:raise ValueError('Replay data hash mismatch')
        records.append((m,s))
    speedmax=max(float(np.linalg.norm(s[1:,4:7],axis=1).max()*1e3) for m,s in records)
    clim=(0,math.ceil(speedmax*2)/2)
    p=pv.Plotter(off_screen=True,window_size=(1600,1000));p.set_background('#090d16')
    p.add_mesh(wall,color='#afc9e5',opacity=.30,smooth_shading=True,specular=.35,show_scalar_bar=False)
    for m,s in records:
        xyz=(s[:,1:4]-origin)*1e6;line=pv.PolyData(xyz);line.lines=np.r_[len(xyz),np.arange(len(xyz))]
        # First row is a free-field diagnostic, not an accepted trajectory speed.
        speeds=np.linalg.norm(s[:,4:7],axis=1)*1e3
        if len(s)>1:speeds[0]=speeds[1]
        line['Track speed (mm/s)']=speeds
        p.add_mesh(line,scalars='Track speed (mm/s)',cmap='turbo',clim=clim,line_width=2.1,opacity=.78,
            scalar_bar_args=dict(title='Track speed (mm/s)',color='#edf1fa',vertical=True,position_x=.88,position_y=.27,height=.48,width=.07,title_font_size=16,label_font_size=14,n_labels=5))
    p.show_grid(grid='back',location='outer',xtitle='X (um)',ytitle='Y (um)',ztitle='Z (um)',
        color='#a7b7ca',font_size=12,n_xlabels=4,n_ylabels=4,n_zlabels=5,all_edges=True,
        fmt='%.0f',show_xlabels=True,show_ylabels=True,show_zlabels=True)
    p.add_text('2.0 mm/s + P9-A | Smoke trajectories',position=(32,948),font_size=22,color='#f1f4fa')
    p.add_text('12 smoke tracks | Production gate failed',position=(32,905),font_size=14,color='#ffb0aa')
    p.add_text('Coordinates relative to vessel bounding-box center (um)',position=(32,8),font_size=10,color='#b5c4d7')
    radius=float(np.linalg.norm(wall.points,axis=1).max());radial=np.linalg.norm(wall.points[:,:2],axis=1)
    elev=np.deg2rad(18);scale=max(np.max(abs(wall.points[:,2])*np.cos(elev)+radial*np.sin(elev)),radial.max()/1.4)/.82
    p.enable_parallel_projection()
    def camera(deg):
        theta=np.deg2rad(deg);direction=np.array([np.cos(elev)*np.cos(theta),np.cos(elev)*np.sin(theta),np.sin(elev)])
        p.camera.position=radius*4*direction;p.camera.focal_point=(0,0,0);p.camera.up=(0,0,1);p.camera.parallel_scale=scale
    camera(35)
    csvrows=[]
    for m,s in records:
        for row in s:csvrows.append([m['particle_id'],*row[:4]])
    with (REPORT/'data/figure04_smoke_xyz.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['particle_id','elapsed_time_s','x_m','y_m','z_m']);w.writerows(csvrows)
    p.screenshot(str(FIG/'04_new_smoke_trajectories.png'))
    movie=ANI/'01_new_flow_smoke_replay.mp4';fps=24;frames=240;end_age=.085
    writer=imageio_ffmpeg.write_frames(str(movie),(1600,1000),fps=fps,codec='libx264',quality=8,macro_block_size=2)
    writer.send(None)
    p.add_text('Age-aligned replay; not simultaneous concentration',position=(32,70),font_size=12,color='#b5c4d7')
    sphere=pv.Sphere(radius=1,theta_resolution=16,phi_resolution=12)
    slerps={m['particle_id']:Slerp(s[:,0],Rotation.from_quat(s[:,[11,12,13,10]])) for m,s in records if len(s)>1}
    for frame in range(frames):
        age=end_age*frame/(frames-1);camera(35+360*frame/frames);centers=[];radii=[];speeds=[];marks=[]
        for m,s in records:
            if age>s[-1,0]:continue
            k=min(int(np.searchsorted(s[:,0],age,side='right'))-1,len(s)-2)
            f=(age-s[k,0])/(s[k+1,0]-s[k,0]);pos=(1-f)*s[k,1:4]+f*s[k+1,1:4];xyz=(pos-origin)*1e6
            centers.append(xyz);radii.append(m['radius_m']*1e6);speeds.append(float(np.linalg.norm(s[k+1,4:7])*1e3))
            axis=slerps[m['particle_id']](age).apply([1,0,0]);marks.extend([xyz-axis*m['radius_m']*1e6*1.03,xyz+axis*m['radius_m']*1e6*1.03])
        if centers:
            cloud=pv.PolyData(np.array(centers));cloud['radius_um']=radii;cloud['Track speed (mm/s)']=speeds
            glyph=cloud.glyph(geom=sphere,scale='radius_um',orient=False)
            p.add_mesh(glyph,name='moving_mb',scalars='Track speed (mm/s)',cmap='turbo',clim=clim,smooth_shading=True,show_scalar_bar=False)
            lines=pv.PolyData(np.array(marks));lines.lines=np.column_stack([np.full(len(marks)//2,2),np.arange(len(marks)).reshape(-1,2)]).ravel()
            p.add_mesh(lines,name='rotation_markers',color='white',line_width=2.2,show_scalar_bar=False)
        p.add_text(f'Trajectory age: {age*1000:5.1f} ms',position=(32,860),name='age',font_size=16,color='#edf1fa')
        writer.send(np.ascontiguousarray(p.screenshot(return_img=True)[:,:,:3]))
    writer.close();p.close()
    dump(REPORT/'data/animation_manifest.json',dict(path=str(movie.relative_to(REPORT)),sha256=sha256(movie),
        fps=fps,frames=frames,duration_s=frames/fps,represented_trajectory_age_s=[0,end_age],
        coordinate_origin_m=origin.tolist(),source_particle_ids=[m['particle_id'] for m,s in records],
        source_sample_sha256={str(m['particle_id']):m['samples_sha256'] for m,s in records},
        replay='AGE_ALIGNED_INDEPENDENT_TRAJECTORIES_NOT_PHYSICAL_CONCENTRATION',
        stop_behavior='HIDE_AFTER_RECORDED_TERMINATION_NO_EXTRAPOLATION',rotation='SAVED_QUATERNION_SLERP_MARKER',
        camera='FIXED_Z_AXIS_ROTATION_Z_UP',dataset_role='SMOKE_FAILED_NOT_FORMAL_PRODUCTION'))


if __name__=='__main__':
    scientific_plots();replay_visuals()
