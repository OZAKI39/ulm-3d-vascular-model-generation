"""Twelve auditable P8.2 figures; no scientific data are invented for display."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from .particle8_replay import read,digest
from .particle81_replay import TITLE,FOOTERS,OUTLETS
from .particle81_visuals import VesselView,COLORS
from .particle82_results import Scene
from .particle82_diagnostics import basin_rows
from .particle82_provenance import atomic_json


def canvas(title,subtitle=''):
    fig=plt.figure(figsize=(16,10),dpi=120,facecolor='white')
    fig.text(.035,.955,'Particle-8.2 | '+title,fontsize=22,weight='bold',color='#183449')
    fig.text(.035,.918,TITLE,fontsize=12,color='#208ca6')
    fig.text(.035,.883,subtitle,fontsize=10,color='#183449')
    for i,line in enumerate(FOOTERS):fig.text(.035,.087-i*.023,line,fontsize=8,color='#884b43')
    return fig


def save(fig,root,name,scene,sources):
    root=Path(root);(root/'figures').mkdir(exist_ok=True);(root/'figure_sources').mkdir(exist_ok=True)
    fig.savefig(root/'figures'/(name+'.png'));plt.close(fig)
    atomic_json(root/'figure_sources'/(name+'.json'),dict(source_scene_sha256=scene.sha256,
        diagnostic_source_sha256={str(p):digest(p) for p in sources},renderer_no_rng=True,
        labels=[TITLE,*FOOTERS],figure=name))


def scale_bar(image,view):
    from PIL import ImageDraw
    from .particle81_visuals import text
    draw=ImageDraw.Draw(image);y=image.height-30
    draw.line((25,y,25+20*image.height/(2*view.scale),y),fill='#183449',width=3)
    text(draw,(25,y-25),'20 µm | equal physical scale',14)
    return image


def generate(root,baseline_points):
    root=Path(root);baseline_points=Path(baseline_points);scene=Scene(root)
    scaling_path=root/'scaling/SERVER_SCALING_BENCHMARK.json';scaling=read(scaling_path)
    admission_path=root/'admission/ADMISSION_BASIN_AUDIT.json';admission=read(admission_path)
    stop_path=root/'admission/stop_audit/SAFETY_STOP_AUDIT.json';stop=read(stop_path)
    extended_path=root/'continuation/EXTENDED_GUARD_AUDIT.json';extended=read(extended_path)
    timestep_path=root/'timestep/TIMESTEP_AUDIT.json';timestep=read(timestep_path)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig=canvas('00  Server CPU scaling','Fixed 64 original natural IDs; all outcomes retained; Linux fork / read-only Frozen arrays')
    rows=scaling['worker_scaling_results'];workers=[r['workers'] for r in rows]
    ax=fig.add_axes([.08,.23,.38,.59]);ax.plot(workers,[r['tracks_per_hour'] for r in rows],'o-',color='#208ca6')
    ax.set(xlabel='Worker processes',ylabel='Integrated trajectories / hour',xticks=workers)
    ax.axvline(scaling['chosen_workers'],color='#c76b46',linestyle='--',label=f"Selected: {scaling['chosen_workers']}");ax.legend()
    ax=fig.add_axes([.57,.23,.36,.59]);ax.plot(workers,[r['peak_process_pss_sum_bytes']/2**30 for r in rows],'o-',label='Process-tree PSS')
    ax.plot(workers,[r['peak_process_rss_sum_bytes']/2**30 for r in rows],'s--',label='RSS sum (shared pages repeated)')
    ax.set(xlabel='Worker processes',ylabel='Peak memory / GiB',xticks=workers);ax.legend()
    fig.text(.08,.15,f"Speedup {scaling['parallel_speedup']:.2f}x | efficiency {scaling['parallel_efficiency']:.1%} | GPU not used by trajectory solver",fontsize=12)
    save(fig,root,'00_server_compute_scaling',scene,[scaling_path])
    fig=canvas('01  Parallel versus serial parity','Every accepted state byte, birth record, terminal reason, outlet and path statistic compared on the server')
    matrix=np.array([[r['exact_scientific_record_and_sample_bytes'] for r in group['ids']] for group in scaling['parallel_vs_serial']])
    ax=fig.add_axes([.09,.30,.82,.46]);ax.imshow(matrix,aspect='auto',vmin=0,vmax=1,cmap=matplotlib.colors.ListedColormap(['#b7443e','#21887a']))
    ax.set(yticks=np.arange(len(matrix)),yticklabels=[f"{p['workers']} vs 1 worker" for p in scaling['parallel_vs_serial']],xlabel='Fixed stable ID (1–64)')
    ax.set_xticks(np.arange(0,64,7),np.arange(1,65,7));fig.text(.10,.19,'Green = exact equality. Worker PID and timing stay in separate execution receipts.',fontsize=13)
    save(fig,root,'01_parallel_vs_serial_parity',scene,[scaling_path])
    point_rows=basin_rows(baseline_points);positions=np.load(baseline_points/'birth_positions.npz')['positions']
    centered=positions-positions.mean(0);_,_,vt=np.linalg.svd(centered,full_matrices=False);uv=centered@vt[:2].T*1e6
    fig=canvas('02  Inlet → outlet point-tracer basins','100,000 zero-radius flux-weighted diagnostic tracers; unresolved paths remain in the denominator')
    ax=fig.add_axes([.08,.21,.48,.62]);labels=np.array([r['outlet'] or 'UNRESOLVED' for r in point_rows])
    for b,color in [*COLORS.items(),('UNRESOLVED','#999999')]:
        take=labels==b;ax.scatter(uv[take,0],uv[take,1],s=1,c=color,label=f'{b}: {take.sum()}',rasterized=True)
    ax.set(xlabel='Inlet-plane display coordinate 1 / µm',ylabel='Inlet-plane display coordinate 2 / µm',aspect='equal');ax.legend(markerscale=5,fontsize=10)
    summary=read(baseline_points/'POINT_TRACER_SUMMARY.json');ax=fig.add_axes([.67,.29,.28,.45]);xx=np.arange(3)
    observed=[summary['fraction_vs_frozen'][o]['observed_fraction'] for o in OUTLETS]
    expected=[summary['fraction_vs_frozen'][o]['frozen_flux_fraction'] for o in OUTLETS]
    bounds=np.array([summary['fraction_vs_frozen'][o]['simultaneous_bonferroni_exact_interval'] for o in OUTLETS])
    ax.errorbar(xx,observed,yerr=[np.array(observed)-bounds[:,0],bounds[:,1]-observed],fmt='o',label='Point fraction + exact CI')
    ax.scatter(xx,expected,marker='x',s=70,color='#c76b46',label='Frozen Q_OUT / Q_IN');ax.set(xticks=xx,xticklabels=['01','02','03'],ylabel='Fraction of all proposals',xlabel='OUTLET');ax.legend(fontsize=9)
    fig.text(.59,.19,'Bonferroni familywise 95% binomial intervals.\nUnresolved point paths are not silently reassigned.',fontsize=11,color='#884b43')
    save(fig,root,'02_inlet_outlet_basin_map',scene,[baseline_points/'POINT_TRACER_SUMMARY.json'])
    basins=[*OUTLETS,'UNRESOLVED_POINT_PATH'];short=['01','02','03','unresolved']
    fig=canvas('03  Finite-size admission by point basin','Every original P8.1 retry retained; original FiniteSizeAdmission re-executed without moving points or resampling sizes')
    for j,(key,title) in enumerate([('proposal_count','All proposals / retries'),('admitted_count','Actually admitted positions'),('guard_exhausted_by_first_proposal_count','Exhausted IDs: first-proposal basin')]):
        ax=fig.add_axes([.07+j*.31,.25,.25,.55]);vals=[admission['by_basin'][b].get(key,0) for b in basins]
        ax.bar(short,vals,color=[*COLORS.values(),'#999999']);ax.set_title(title,fontsize=11);ax.set_ylabel('Count')
        for i,v in enumerate(vals):ax.text(i,v,str(v),ha='center',va='bottom',fontsize=10)
    fig.text(.07,.17,admission['accounting'][:165]+'…',fontsize=9)
    save(fig,root,'03_admission_bias_by_outlet_basin',scene,[admission_path,root/'admission/ORIGINAL_ADMISSION_RECHECK.json'])
    fig=canvas('04  Diameter conditioning by point basin','Scheduled/exhausted IDs use first proposals; admitted sizes use accepted-position basins; no diameter redraw')
    for j,b in enumerate(basins):
        ax=fig.add_axes([.065+j*.235,.24,.20,.57]);data=admission['diameters'][b]
        for key,label,color in [('scheduled_first_proposal_diameter_um','scheduled','#999999'),('admitted_diameter_um','admitted','#208ca6'),('guard_exhausted_diameter_um','exhausted','#c76b46')]:
            vals=data.get(key,[])
            if vals:ax.hist(vals,bins=np.linspace(.75,5.25,25),histtype='step',linewidth=1.7,label=f'{label} n={len(vals)}',color=color)
        ax.set(title=short[j],xlabel='Original diameter / µm',ylabel='IDs');ax.legend(fontsize=8)
        if not data.get('admitted_diameter_um'):ax.text(.08,.8,'NO ADMITTED CASES',transform=ax.transAxes,color='#884b43',fontsize=9)
    save(fig,root,'04_size_bias_by_outlet_basin',scene,[admission_path])
    # Remaining figures use the same final scene plus their named diagnostics.
    other_figures(root,scene,stop,extended,timestep,[stop_path,extended_path,timestep_path])


def other_figures(root,scene,stop,extended,timestep,paths):
    import pyvista as pv
    fig=canvas('05  Original 720 computational stop locations','Full unchanged Frozen vessel; marker color = point-tracer basin, not a claimed finite-size outlet')
    view=VesselView(scene,size=(1050,650),wall_opacity=.13)
    rows=stop['rows'];xyz=np.array([[r['x_m'],r['y_m'],r['z_m']] for r in rows]);mesh=pv.PolyData((xyz-view.origin)*1e6)
    mesh['basin']=np.array([OUTLETS.index(r['point_tracer_basin'])+1 if r['point_tracer_basin'] in OUTLETS else 0 for r in rows])
    view.plotter.add_mesh(mesh,scalars='basin',cmap=['#999999',*COLORS.values()],clim=(0,3),point_size=5,show_scalar_bar=False,render=False)
    image,_=view.capture(70);scale_bar(image,view);view.close();ax=fig.add_axes([.035,.16,.66,.69]);ax.imshow(image);ax.axis('off')
    ax=fig.add_axes([.75,.32,.20,.41]);counts=stop['by_basin'];ax.bar(['01','02','03','?'],[counts.get(b,0) for b in OUTLETS+['UNRESOLVED_POINT_PATH']],color=[*COLORS.values(),'#999999']);ax.set(ylabel='Original safety stops',xlabel='Point basin')
    fig.text(.73,.21,'Markers are display symbols.\nNo physiological capture claim.\nRaw gaps, times and provider counts saved.',fontsize=11,color='#884b43')
    save(fig,root,'05_safety_stop_spatial_map',scene,[paths[0]])
    fig=canvas('06  Extended computational-guard outcomes','Same initial states, radii, dt and P6.5 acceptance rule; exact original accepted-prefix comparisons')
    for j,trial in enumerate(extended['factors']):
        ax=fig.add_axes([.07+j*.46,.30,.39,.49]);basins=OUTLETS+['UNRESOLVED_POINT_PATH'];x=np.arange(4)
        complete=[sum(r['completed'] and r['point_basin']==b for r in trial['rows']) for b in basins]
        stopped=[sum(r['end_reason']=='INTEGRATION_SAFETY_STOP' and r['point_basin']==b for r in trial['rows']) for b in basins]
        other=[sum(not r['completed'] and r['end_reason']!='INTEGRATION_SAFETY_STOP' and r['point_basin']==b for r in trial['rows']) for b in basins]
        ax.bar(x,complete,label='Stop → complete',color='#208ca6')
        ax.bar(x,stopped,bottom=complete,label='Still safety-stopped',color='#999999')
        ax.bar(x,other,bottom=np.array(complete)+stopped,label='Other terminal reason',color='#c76b46')
        ax.set(xticks=x,xticklabels=['01','02','03','?'],xlabel='Original point basin',ylabel='Original safety-stop IDs',title=f'Guard ×{trial["guard_factor"]}')
        ax.legend(fontsize=8)
        for i,values in enumerate(zip(complete,stopped,other)):ax.text(i,sum(values)+7,' / '.join(map(str,values)),ha='center',fontsize=9)
    fig.text(.075,.19,'Baseline completed: 1249 (1248 from point basin 02, 1 from basin 01). Original stopped: 720.\nOriginal accepted prefixes remain exact; new suffixes start at saved states. No physiological capture claim.',fontsize=11,color='#183449')
    save(fig,root,'06_extended_guard_outcomes',scene,[paths[1]])
    fig=canvas('07  Basin-stratified timestep audit','dt / dt/2 / dt/4 with guard ×16; NOT PRODUCTION TIMESTEP SELECTION')
    for j,basin in enumerate(OUTLETS):
        ax=fig.add_axes([.075+j*.305,.27,.25,.53]);group=timestep['groups'][basin]
        rows=[r for r in timestep['rows'] if r['point_basin']==basin]
        if rows:
            for div in [1,2,4]:
                rr=[r for r in rows if r['dt_divisor']==div];ax.scatter(np.full(len(rr),div),[r['trajectory_deviation_max_m']*1e6 for r in rr],s=15,alpha=.7,color=COLORS[basin])
            ax.set(xlabel='dt divisor',ylabel='Max common-time path deviation / µm',xticks=[1,2,4])
        else:ax.text(.08,.55,'NO ADMITTED CASES\nNo artificial basin replacements',transform=ax.transAxes,color='#884b43');ax.set_xticks([]);ax.set_yticks([])
        ax.set_title(f'{basin}: available {group["available"]}',fontsize=11)
    save(fig,root,'07_multioutlet_timestep_audit',scene,[paths[2]])
    coverage=scene.catalog['all_outlets_observed'];name='08_full_network_ulm_reconstruction' if coverage else '08_observed_network_coverage'
    fig=canvas('08  '+('Full-network' if coverage else 'Observed-network')+' ULM-style trajectories','Complete natural MB paths only; unobserved branches remain grey; no acoustic image-formation model')
    view=VesselView(scene,size=(720,600));view.paths(scene.completed,'ensemble',opacity=.6,width=1.2);image,_=view.capture(65);scale_bar(image,view);view.close()
    ax=fig.add_axes([.025,.20,.44,.64]);ax.imshow(image);ax.axis('off')
    wall=np.load(scene.root/'data/full_frozen_geometry.npz')['WALL_points_m']*1e6
    for i,(axes,label) in enumerate([((0,1),'XY'),((0,2),'XZ'),((1,2),'YZ')]):
        ax=fig.add_axes([.52+(i%2)*.23,.55-(i//2)*.35,.20,.26]);ax.scatter(wall[:,axes[0]],wall[:,axes[1]],s=.1,c='#dddddd')
        for o in OUTLETS:
            segments=[scene.arrays[pid][:,1:4][:,axes]*1e6 for pid in scene.completed if scene.entries[pid]['exit_outlet']==o]
            if segments:ax.add_collection(LineCollection(segments,colors=COLORS[o],linewidths=.25,alpha=.3))
        ax.autoscale();ax.set_aspect('equal');ax.set(title=label,xlabel='µm',ylabel='µm');ax.tick_params(labelsize=8)
    fig.text(.77,.25,'Natural completed paths\n'+ '\n'.join(f'{o}: {scene.catalog["outlet_counts"][o]}' for o in OUTLETS),fontsize=12)
    save(fig,root,name,scene,[])
    fig=canvas('09  Natural finite-size MB outlet statistics','Natural flux weighting preserved; failures retained; outlet proportions are not uncensored physiological probabilities')
    ax=fig.add_axes([.08,.26,.39,.54]);values=[scene.catalog['outlet_counts'][o] for o in OUTLETS];ax.bar(['01','02','03'],values,color=list(COLORS.values()));ax.set(xlabel='OUTLET',ylabel='Completed natural MBs')
    for i,v in enumerate(values):ax.text(i,v,str(v),ha='center',va='bottom')
    ax=fig.add_axes([.57,.26,.36,.54]);c=scene.catalog;vals=[c['completed'],c['admitted']-c['completed'],c['scheduled']-c['admitted']]
    ax.bar(['complete','stopped','unadmitted'],vals,color=['#208ca6','#999999','#c76b46']);ax.set_ylabel('All scheduled natural IDs')
    for i,v in enumerate(vals):ax.text(i,v,str(v),ha='center',va='bottom')
    save(fig,root,'09_outlet_statistics',scene,[])
    fig=canvas('11  Scientific and compute limitations remain visible')
    lines=[('FROZEN, INDEPENDENT TRAJECTORIES','No RBC-resolved suspension, real RBC passage, MB–MB coupling or complete PK.'),
        ('POINT-TRACER CENSORING','Unresolved point paths remain explicit; a boundary-flux match is not forced.'),
        ('FINITE-SIZE ADMISSION','Accepted MB locations/sizes are conditioned by the unchanged original wall rule.'),
        ('COMPUTATIONAL STOPS','Extra guard budget cannot prove physiological capture or downstream passage.'),
        ('NETWORK COVERAGE','Only observed natural complete tracks reconstruct branches; no copied/synthetic paths.'),
        ('NUMERICAL SENSITIVITY','Production timestep, neighbors and non-spherical lubrication remain NOT FROZEN.'),
        ('COMPUTE PROVENANCE','P8.1 remote heavy compute NOT_PROVEN; P8.2 formal receipts checked individually.'),
        ('DISPLAY AND REVIEW','Camera/interpolation/time cuts are display-only. User manual review remains pending.')]
    for i,(title,detail) in enumerate(lines):
        y=.81-i*.082;fig.text(.075,y,title,fontsize=13,color='#208ca6',weight='bold');fig.text(.075,y-.03,detail,fontsize=11)
    save(fig,root,'11_particle8_2_limitations',scene,[])
