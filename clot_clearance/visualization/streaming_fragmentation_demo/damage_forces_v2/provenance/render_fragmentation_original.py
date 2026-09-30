"""Fixed-camera, true-position topology/transport visualization from saved data."""
import argparse,json,os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import hsv_to_rgb,Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import imageio.v2 as imageio
import imageio_ffmpeg


def shown_edges(pairs,active,X):
    """Spanning forest guarantees displayed connectivity, plus short local bonds."""
    n=len(X);parent=np.arange(n);keep=[]
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    candidate=np.flatnonzero(active)
    order=candidate[np.argsort(np.linalg.norm(X[pairs[candidate,1]]-X[pairs[candidate,0]],axis=1),kind='stable')]
    for k in order:
        a,b=pairs[k];a=find(a);b=find(b)
        if a!=b:parent[a]=b;keep.append(k)
    return np.asarray(keep,dtype=int)


def color_ids(ids):
    hue=np.mod(ids*.61803398875+.06,1)
    return hsv_to_rgb(np.column_stack((hue,np.full(len(ids),.65),np.full(len(ids),.96))))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();run=args.run.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'frames').mkdir();(out/'figures').mkdir()
    z=np.load(run/'states.npz');h=json.loads((run/'history.json').read_text());summary=json.loads((run/'SUMMARY.json').read_text())
    config=json.loads((run/'CONFIG.json').read_text());comps=json.loads((run/'components.json').read_text())
    x=z['x']*1e3;X=z['X']*1e3;pairs=z['pairs'];fixed=z['fixed'];bubble=np.array(config['streaming']['bubble_center_m'])*1e3
    xclear=config['transport']['x_clearance_m']*1e3
    xlim=(min(X[:,0].min()-.15,-.8),max(x[:,:,0].max()+.15,xclear+.25));ylim=(-.55,.55)
    zlim=(min(x[:,:,2].min()-.12,-1.2),max(x[:,:,2].max()+.15,bubble[2]+.2))
    dark='#07101e';white='#e8f0fa';muted='#aabbd0';blue='#55c7ff'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    forests=[shown_edges(pairs,b,X) for b in z['bond_active']]

    def scene(ax,k,mode,dark_mode=True,annotate=True):
        fg=white if dark_mode else '#1f2937';ax.set_facecolor(dark if dark_mode else 'white')
        active=forests[k];segments=x[k,pairs[active]]
        ax.add_collection3d(Line3DCollection(segments,colors='#69a7c5' if dark_mode else '#6d8799',linewidths=.45,alpha=.5))
        recent=np.flatnonzero(z['recent_broken'][k])
        # For dense bursts retain a deterministic sample for readable fracture markers.
        if len(recent)>350:recent=recent[np.linspace(0,len(recent)-1,350,dtype=int)]
        if len(recent):ax.add_collection3d(Line3DCollection(x[k,pairs[recent]],colors='#fb5d65',linewidths=.65,alpha=.6))
        colors=colormaps['inferno'](z['damage'][k]) if mode=='damage' else color_ids(z['fragment_id'][k])
        mobile=~fixed
        ax.scatter(*x[k,mobile].T,c=colors[mobile],s=13 if dark_mode else 10,edgecolors='#9daec1' if dark_mode else '#374151',linewidths=.2,depthshade=False,zorder=5)
        ax.scatter(*x[k,fixed].T,c='#9fa9b5',marker='s',s=14,depthshade=False,zorder=4)
        ax.scatter(*bubble,c=blue,s=70,edgecolors=fg,linewidths=.5,depthshade=False,zorder=15)
        # Marker is enlarged only; all particle coordinates use true displacement.
        ax.plot([xclear,xclear],[ylim[0],ylim[0]],[zlim[0],zlim[1]],'--',color='#55df9b',lw=1.2)
        ax.plot([xclear,xclear],[ylim[0],ylim[1]],[zlim[0],zlim[0]],'--',color='#55df9b',lw=1.2)
        ax.quiver(xlim[0]+.08,ylim[0],zlim[1]-.05,.5,0,0,color=blue,arrow_length_ratio=.18,linewidth=1.7)
        ax.text(xlim[0]+.05,ylim[0],zlim[1]+.025,'Fluid +X',color=blue,fontsize=9)
        ax.set(xlim=xlim,ylim=ylim,zlim=zlim,xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)')
        ax.set_box_aspect((xlim[1]-xlim[0],ylim[1]-ylim[0],zlim[1]-zlim[0]));ax.view_init(elev=22,azim=-68)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.label.set_color(fg);axis.line.set_color('#5c7086' if dark_mode else '#9aa7b5')
            axis.set_pane_color((.04,.075,.12,1) if dark_mode else (.97,.98,.99,1))
            axis._axinfo['grid']['color']='#293c51' if dark_mode else '#dde3e9'
        ax.tick_params(colors=fg,labelsize=8,pad=1)
        if annotate and mode=='fragments':
            mobile_components=sorted([a for a in comps[k]['components'] if not a['attached_to_base']],key=lambda a:-a['particle_count'])[:1]
            for comp in mobile_components:
                pos=np.array(comp['center_of_mass_m'])*1e3
                ax.text(pos[0],pos[1],pos[2]+.08,f'ID {comp["fragment_id"]}',color=fg,fontsize=8,zorder=20)

    events=summary['events'];B=events['first_bond_failure']['macro_step'];C=events['first_detachment']['macro_step'];D=len(h)-1
    selected=[0,B,C,D]
    fig=plt.figure(figsize=(16,10),layout='constrained')
    for j,(k,title) in enumerate(zip(selected,['A  Intact clot','B  First bond failure','C  First detachment','D  Fragment transport'])):
        ax=fig.add_subplot(2,2,j+1,projection='3d',computed_zorder=False);scene(ax,k,'damage' if j<2 else 'fragments',False)
        row=h[k]
        ax.set_title(f'{title} · N = {row["represented_cycles"]:,}\nMean/max damage {row["mean_damage"]:.3f}/{row["maximum_damage"]:.3f} · broken {100*row["broken_bond_fraction"]:.2f}%\nAttached {100*row["attached_volume_fraction"]:.1f}% · detached {100*row["detached_volume_fraction"]:.1f}% · free components {row["detached_components"]}',fontsize=11,pad=10)
    fig.suptitle('Load-driven PD failure → detachment → prescribed-fluid transport\nTrue positions; active spanning-forest bonds in blue-gray, newly broken bonds in red',fontsize=16)
    fig.savefig(out/'figures/four_panel_summary.png',dpi=300);fig.savefig(out/'figures/four_panel_summary.pdf');plt.close(fig)
    keys=['mean_damage','maximum_damage','broken_bond_fraction','attached_volume_fraction','detached_volume_fraction','cleared_volume_fraction','detached_components']
    labels=['Mean particle damage','Maximum particle damage','Broken bond fraction','Attached volume fraction','Detached volume fraction','Cleared volume fraction','Detached component count']
    fig,axes=plt.subplots(2,4,figsize=(16,8),layout='constrained');N=z['cycles']
    event_names=['first_bond_failure','first_detachment','first_fragment','first_clearance_crossing']
    event_labels=['First bond failure','First detachment','First fragment','First clearance crossing'];event_colors=['#cf3b50','#7551ae','#c79522','#229666']
    for ax,key,label in zip(axes.ravel(),keys,labels):
        ax.plot(N,[row[key] for row in h],'-o',lw=1.6,ms=3,color='#1879ae');ax.set(xlabel='Represented cycles N',ylabel=label);ax.grid(alpha=.22)
        for e,color,linestyle in zip(event_names,event_colors,['--','--',':','--']):
            if events[e]:ax.axvline(events[e]['cycles'],color=color,lw=1.2,ls=linestyle)
    legend=axes.ravel()[-1];legend.axis('off')
    for i,(e,label,color) in enumerate(zip(event_names,event_labels,event_colors)):
        legend.text(0,.9-i*.17,f'{label}\nN = {events[e]["cycles"]:,}' if events[e] else f'{label}: none',color=color,fontsize=12)
    legend.text(0,.05,'Qualitative forced-failure verification\nUncalibrated damage and transport',fontsize=11,color='#49576b')
    fig.suptitle('Damage, connectivity and clearance histories',fontsize=17)
    fig.savefig(out/'figures/time_histories.png',dpi=300);fig.savefig(out/'figures/time_histories.pdf');plt.close(fig)
    os.environ['IMAGEIO_FFMPEG_EXE']=imageio_ffmpeg.get_ffmpeg_exe()
    for focus in ['damage','fragments']:
        folder=out/'frames'/focus;folder.mkdir();images=[]
        for k,row in enumerate(h):
            fig=plt.figure(figsize=(16,9),facecolor=dark)
            modes=[focus,'fragments' if focus=='damage' else 'damage']
            for j,mode in enumerate(modes):
                ax=fig.add_axes([.015+.49*j,.24,.48,.60],projection='3d',facecolor=dark,computed_zorder=False)
                scene(ax,k,mode,True)
                ax.set_title('Particle damage' if mode=='damage' else 'Stable fragment IDs',color=white,fontsize=15,pad=3)
            fig.text(.035,.955,'STREAMING → FAILURE → FRAGMENT TRANSPORT',color=white,fontsize=24,weight='bold')
            fig.text(.035,.900,f'N = {row["represented_cycles"]:,} represented cycles   |   fixed camera · true displacement scale',color=blue,fontsize=17)
            fig.text(.04,.17,f'Mean/max damage  {row["mean_damage"]:.3f} / {row["maximum_damage"]:.3f}\nBroken bonds         {100*row["broken_bond_fraction"]:.2f}%',color=white,fontsize=14,va='top',linespacing=1.6)
            fig.text(.38,.17,f'Attached / detached   {100*row["attached_volume_fraction"]:.1f}% / {100*row["detached_volume_fraction"]:.1f}%\nDetached components   {row["detached_components"]}',color=white,fontsize=14,va='top',linespacing=1.6)
            fig.text(.73,.17,f'Cleared volume   {100*row["cleared_volume_fraction"]:.1f}%\nMaximum travel   {row["maximum_fragment_travel_m"]*1e3:.3f} mm',color=white,fontsize=14,va='top',linespacing=1.6)
            fig.text(.04,.055,'Gray squares: fixed base · blue marker: prescribed bubble (enlarged) · green line: clearance plane · red: newly broken bonds',color=muted,fontsize=11)
            fig.text(.04,.025,'All particles retained. Synthetic loading and empirical failure/transport parameters; not a quantitative sonothrombolysis prediction.',color=muted,fontsize=11)
            cax=fig.add_axes([.40,.225,.20,.012]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='inferno'),cax=cax,orientation='horizontal')
            cb.ax.tick_params(colors=white,labelsize=8,pad=1);cb.set_label('Particle damage',color=white,fontsize=9,labelpad=1)
            path=folder/f'frame_{k:04d}.png';fig.savefig(path,dpi=100,facecolor=dark);plt.close(fig);images.append(imageio.imread(path)[:,:,:3])
        with imageio.get_writer(out/f'fragmentation_{focus}.mp4',fps=12,codec='libx264',quality=8,macro_block_size=1) as writer:
            for im in images:
                for _ in range(6):writer.append_data(im)
        imageio.mimsave(out/f'fragmentation_{focus}.gif',images,duration=500,loop=0)
    info=dict(run=str(run),frame_count=len(h),fps=12,repeats_per_saved_state=6,movie_seconds=len(h)/2,
        selected_event_frame_indices=selected,camera=dict(elevation_deg=22,azimuth_deg=-68,xlim_mm=xlim,ylim_mm=ylim,zlim_mm=zlim),
        displacement_scale=1,equal_axis_units=True,active_bond_display='spanning forest of ALL active graph components',recent_broken_display='deterministic sample up to 350 per macro',
        particle_display='all particles; no culling or deletion',damage_color_range=[0,1])
    (out/'RENDER_MANIFEST.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info,indent=2))


if __name__=='__main__':main()
