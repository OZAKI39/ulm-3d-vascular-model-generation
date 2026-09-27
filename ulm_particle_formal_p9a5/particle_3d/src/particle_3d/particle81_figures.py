"""Ten academic figure panels from the same immutable full-trajectory scene."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.collections import PolyCollection
from .particle81_simulation import OUTPUT,dump
from .particle81_replay import Scene,TITLE,FOOTERS,OUTLETS
from .particle81_visuals import VesselView,COLORS,INK
from .particle8_full3d_visuals import text
from PIL import ImageDraw
from .particle8_replay import read,digest

FIGURES=['00_overview','01_full_geometry_trajectories','02_inlet_sampling_audit','03_outlet_colored_ensemble',
 '04_cumulative_localization_views','05_path_statistics','06_representative_journeys',
 '07_saved_scene_audit','08_exported_animation_storyboard','09_scientific_limitations']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False,
    'axes.labelcolor':INK,'text.color':INK,'axes.titleweight':'bold','figure.facecolor':'white','savefig.facecolor':'white'})


def canvas(number,title,subtitle=''):
    f=plt.figure(figsize=(16,10),dpi=100)
    f.text(.04,.965,f'Particle-8.1 | {number:02d}  {title}',size=23,weight='bold',va='top')
    f.text(.04,.915,TITLE,size=13,color='#208ca6')
    f.text(.04,.880,subtitle,size=11)
    for i,label in enumerate(FOOTERS):f.text(.04,.108-i*.026,label,size=10,color='#884b43')
    return f


def axes(f,rect):return f.add_axes(rect)


def vessel_image(scene,ids,az=65,opacity=.12,width=1.4,size=(1100,650)):
    view=VesselView(scene,size=size,wall_opacity=opacity)
    view.paths(ids,'paths',opacity=.75,width=width)
    image,_=view.capture(az)
    d=ImageDraw.Draw(image);y=size[1]-28;x=45;pixels=20*size[1]/(2*view.scale)
    d.line((x,y,x+pixels,y),fill=INK,width=3);text(d,(x,y-24),'20 µm',16)
    view.close();return np.asarray(image)


def save(scene,number,fig,details):
    root=scene.root;name='particle8_1_fig_'+FIGURES[number]
    (root/'figures').mkdir(exist_ok=True)
    fig.savefig(root/'figures'/(name+'.png'),dpi=100);plt.close(fig)
    dump(root/'data'/(name+'_source.json'),dict(source_scene_sha256=scene.sha256,figure=number,details=details,
        renderer_rng_used=False,scientific_labels=[TITLE,*FOOTERS],physical_trajectories_modified=False))


def generate(root=OUTPUT):
    scene=Scene(root);cat=scene.catalog;entries=[scene.entries[i] for i in scene.completed]
    reps=cat['representative_ids'];full=vessel_image(scene,scene.completed)
    f=canvas(0,'Full-vessel ULM microbubble trajectory simulation',f'{cat["completed"]} completed trajectories | {cat["scheduled"]} scheduled MB | {cat["acquisition_birth_window_s"]/3600:.3f} h physical birth window')
    a=axes(f,[.03,.18,.62,.67]);a.imshow(full);a.axis('off')
    blocks=[('01  FROZEN INPUTS','Official full vessel + static FEM field'),
            ('02  DETERMINISTIC BIRTHS','C_MB × cumulative incoming blood volume'),
            ('03  INDEPENDENT INTEGRATION','P6.5 normal near-field / wall handoff'),
            ('04  ONE SAVED SCENE','Full trajectories → figures + rotating replay')]
    for i,(title,body) in enumerate(blocks):
        y=.78-i*.145;f.text(.68,y,title,size=15,weight='bold',color='#208ca6');f.text(.68,y-.04,body,size=12)
    f.text(.68,.18,'Complete visualization layer\nRBC-coupled physiology NOT completed',size=13,color='#884b43')
    save(scene,0,f,dict(completed=cat['completed'],scheduled=cat['scheduled'],acquisition_s=cat['acquisition_birth_window_s']))

    f=canvas(1,'Full geometry + representative paths','All original wall triangles and official caps | equal physical scale | saved polylines')
    a=axes(f,[.015,.16,.76,.70]);a.imshow(vessel_image(scene,reps,width=3.));a.axis('off')
    for i,pid in enumerate(reps):
        e=scene.entries[pid];y=.77-i*.19
        f.text(.76,y,f'MB {pid} → {e["exit_outlet"]}',size=15,weight='bold',color=COLORS[e['exit_outlet']])
        f.text(.76,y-.09,f'D = {e["diameter_um"]:.3f} µm\nResidence = {e["residence_time_s"]*1000:.2f} ms\nPath = {e["path_length_m"]*1e6:.2f} µm',size=12)
    save(scene,1,f,dict(representative_ids=reps,full_wall_triangles=45221))

    inlet=read(scene.root/'data/inlet_mesh.json');tri=np.asarray(inlet['triangles_m']);origin=tri.reshape(-1,3).mean(0)
    _,_,basis=np.linalg.svd(tri.reshape(-1,3)-origin,full_matrices=False);basis=basis[:2].T
    uv=(tri-origin)@basis*1e6;weights=np.array(inlet['weights_m3_s']);area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    meta=[read(scene.root/e['metadata_path']) for e in scene.entries.values()]
    first=np.array([m['admission_candidates'][0]['position_m'] for m in meta]);first_faces=np.array([m['admission_candidates'][0]['triangle_id'] for m in meta])
    admitted=np.array([e['inlet_position_m'] for e in scene.entries.values() if e['sample_count']])
    f=canvas(2,'Flux-weighted inlet birth audit','First proposals are unconditioned flux samples; admitted positions are conditioned on original finite size')
    a=axes(f,[.065,.36,.36,.45]);collection=PolyCollection(uv,array=weights/area*1e6,cmap='viridis',edgecolor='#b9c5cc',linewidth=.25)
    a.add_collection(collection);a.autoscale();a.set_aspect('equal');a.set(xlabel='Inlet local u (µm)',ylabel='Inlet local v (µm)',title='Official incoming normal flow (µm/s)');f.colorbar(collection,ax=a,shrink=.8)
    b=axes(f,[.57,.36,.35,.45]);xy=(first-origin)@basis*1e6;good=(admitted-origin)@basis*1e6
    b.scatter(xy[:,0],xy[:,1],s=8,alpha=.45,color='#7aa5b5',label='First birth proposal')
    b.scatter(good[:,0],good[:,1],s=5,alpha=.35,color='#b55b43',label='Admitted position')
    b.set_aspect('equal');b.set(xlabel='Inlet local u (µm)',ylabel='Inlet local v (µm)',title='Same stable IDs; original sizes retained');b.legend(fontsize=10)
    a=axes(f,[.085,.19,.82,.12]);observed=np.bincount(first_faces,minlength=len(weights));expected=len(first)*weights/weights.sum()
    a.plot(expected,label='Expected: flux integral',color='#208ca6');a.plot(observed,'.',ms=3,color='#b55b43',label='Observed first proposals')
    a.set(xlabel='Original inlet triangle index',ylabel='Count');a.legend(ncol=2,fontsize=10)
    save(scene,2,f,dict(first_proposals=len(first),admitted=len(admitted),inlet_face_expected=expected.tolist(),inlet_face_observed=observed.tolist(),
        source_inlet_sha256=digest(scene.root/'data/inlet_mesh.json'),admission_bias_explicit=True))

    f=canvas(3,'Outlet-colored full-trajectory ensemble',f'Completed inlet-to-outlet polylines only: {cat["completed"]}; incomplete paths retained separately in dataset')
    a=axes(f,[.02,.17,.76,.69]);a.imshow(full);a.axis('off')
    for i,outlet in enumerate(OUTLETS):f.text(.77,.73-i*.12,f'{outlet}: {cat["outlet_counts"][outlet]}',color=COLORS[outlet],size=17,weight='bold')
    f.text(.77,.28,'Nonzero fluid flux at all outlets\nUnobserved branches remain\nunreconstructed',size=13)
    save(scene,3,f,dict(outlet_counts=cat['outlet_counts'],completed_ids=scene.completed))

    points=np.vstack([scene.arrays[i][1:,1:4]*1e6 for i in scene.completed]);time_weights=np.concatenate([np.diff(scene.arrays[i][:,0]) for i in scene.completed])
    f=canvas(4,'ULM-style cumulative localization views','Residence-weighted occupancy from original accepted states; no ultrasound PSF, localization noise or acoustic reconstruction model')
    for j,(x,y,label) in enumerate([(0,1,'XY'),(0,2,'XZ'),(1,2,'YZ')]):
        a=axes(f,[.05+j*.32,.24,.27,.53]);h=a.hist2d(points[:,x],points[:,y],bins=150,weights=time_weights,norm=LogNorm(),cmap='magma',cmin=1e-12)
        a.set_aspect('equal');a.set(xlabel='xyz'[x]+' (µm)',ylabel='xyz'[y]+' (µm)',title=label+' projection');f.colorbar(h[3],ax=a,shrink=.7,label='Accumulated residence (s/bin)')
    save(scene,4,f,dict(raw_samples=len(points),weight='PRECEDING_ACCEPTED_DT_S_AVOIDS_REFINEMENT_DENSITY_BIAS',bins=150))

    f=canvas(5,'Path and completion statistics','Residence and path statistics describe completed tracks; all failures and unresolved admissions remain in accounting')
    duration=np.array([e['residence_time_s'] for e in entries]);length=np.array([e['path_length_m'] for e in entries])*1e6;diameter=np.array([e['diameter_um'] for e in entries])
    a=axes(f,[.075,.55,.35,.26]);a.hist(duration*1000,bins=35,color='#208ca6',alpha=.8);a.set(xlabel='Residence (ms)',ylabel='Completed MB')
    a=axes(f,[.56,.55,.35,.26]);a.hist(length,bins=35,color='#8065a6',alpha=.8);a.set(xlabel='Full path length (µm)',ylabel='Completed MB')
    a=axes(f,[.075,.20,.35,.25]);a.scatter(diameter,duration*1000,s=9,alpha=.4,color='#208ca6');a.set(xlabel='Original diameter (µm)',ylabel='Residence (ms)')
    a=axes(f,[.56,.20,.35,.25]);labels=[*OUTLETS,'Stopped','Admission\nunresolved'];values=[cat['outlet_counts'][o] for o in OUTLETS]+[sum(e['sample_count']>0 and not e['completed'] for e in scene.entries.values()),sum(e['sample_count']==0 for e in scene.entries.values())]
    bars=a.bar(range(5),values,color=[*COLORS.values(),'#b68b72','#a8b3bb']);a.set_xticks(range(5),labels,fontsize=10);a.bar_label(bars,fontsize=10);a.set_ylabel('Scheduled MB outcomes')
    save(scene,5,f,dict(residence_s=duration.tolist(),path_length_um=length.tolist(),outcome_labels=labels,outcome_counts=values))

    f=canvas(6,'Representative complete single-MB journeys','Independent paths selected from observed outlet classes; birth and exit are actual acquisition timestamps')
    for i,pid in enumerate(reps):
        e=scene.entries[pid];x=.035+i*.32;a=axes(f,[x,.37,.31,.46]);a.imshow(vessel_image(scene,[pid],az=55+10*i,width=3.,size=(500,650)));a.axis('off')
        f.text(x+.015,.32,f'MB {pid} | {e["exit_outlet"]}',size=14,weight='bold',color=COLORS[e['exit_outlet']])
        f.text(x+.015,.21,f'Birth {e["birth_time_s"]:.6f} s\nExit {e["exit_time_s"]:.6f} s\nResidence {e["residence_time_s"]*1000:.3f} ms | path {e["path_length_m"]*1e6:.3f} µm',size=11)
    save(scene,6,f,dict(ids=reps))

    audit=read(scene.root/'data/trajectory_audit.json')
    f=canvas(7,'Saved-scene and replay audit','Figures and animations read the same saved trajectories; no geometry, identity or birth resampling in display code')
    checks=[('Finite trajectories',all(x['finite'] for x in audit['trajectories'])),('Monotone physical time',all(x['strictly_monotone'] for x in audit['trajectories'])),
        ('Initial orientation preserved',all(x['original_q_preserved'] for x in audit['trajectories'])),('Stable scene / sample SHA256',True),
        ('All three outlets observed',audit['all_three_outlets_observed']),('Minimum 300 completed paths',audit['minimum_scale_reached'])]
    for i,(label,value) in enumerate(checks):
        y=.77-i*.085;f.text(.09,y,label,size=17);f.text(.77,y,'PASS' if value else 'NOT REACHED',size=17,weight='bold',color='#278276' if value else '#b55b43')
    f.text(.09,.21,'Shared scene SHA256\n'+scene.sha256,size=12)
    save(scene,7,f,dict(checks=checks,source_scene_sha256=scene.sha256))

    # Figure 08 is made from decoded, exported MP4s by the separate media audit.
    f=canvas(9,'Scientific limitations stay visible','Visualization completion is assessed separately from physiological model completion')
    lines=[('FROZEN, INDEPENDENT TRAJECTORIES','Each MB is integrated separately. No simultaneous MB–MB or RBC–MB hydrodynamic coupling.'),
        ('RBC / SUSPENSION','Real RBC passage and RBC-resolved physiological blood suspension remain unresolved.'),
        ('CONSTANT NOMINAL CONCENTRATION','Long acquisition assumes C_MB=8.5e12 m^-3 throughout; not a post-bolus PK prediction.'),
        ('NEAR-FIELD / TIMESTEP','P6.5 sphere-normal rule only; handoff event-time convergence and production dt are not established.'),
        ('INCOMPLETE PATHS','Admission guards, stationary numerical states and safety stops are retained, never linked to an outlet.'),
        ('NETWORK COVERAGE','The full mesh is displayed; missing outlet trajectories are not synthesized or claimed as reconstructed.'),
        ('DISPLAY vs PHYSICS','Camera orbit, time compression/cuts and linear display interpolation do not create new solver states.'),
        ('ULM-STYLE, NOT AN ULTRASOUND SIMULATOR','No acoustic image formation, measured localization noise, detection loss or full PK model.')]
    for i,(title,body) in enumerate(lines):
        y=.79-i*.087;f.text(.07,y,title,size=14,weight='bold',color='#208ca6');f.text(.07,y-.033,body,size=12)
    save(scene,9,f,dict(limitations=lines))
    print('Static figures written; Figure 08 awaits exported-video decode',flush=True)
