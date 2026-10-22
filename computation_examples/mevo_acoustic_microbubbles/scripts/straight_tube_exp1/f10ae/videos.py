"""True-radius spheres at exact accepted frames; no interpolation or reversal."""
import time,numpy as np,av
from PIL import Image,ImageDraw
from scripts.straight_tube_exp1.f10aa.videos import scene,font
from scripts.straight_tube_exp1.f10ab.metrics import measure
from .figures import traces
from .common import *
def main():
    out=BASE/'videos';out.mkdir(exist_ok=True);rows=read(BASE/'ALL_RUNS.json');by={r['id']:r for r in rows};sel=read(BASE/'SELECTED_FROZEN.json')['selected'];validation=read(BASE/'VALIDATION.json');admission=read(BASE/'INTERPRETATION.json');t1=by[sel['T1']['source']];t4=by[sel['T4']['source']]
    stress=read(BASE/'ROBUSTNESS_RUNS.json');cc=read(BASE/'ROBUSTNESS_CONVERGENCE.json');freq=[]
    for f in [2400000,2500000,2600000]:
        r=next(r for r in stress if r['id']==f'freq_T1_{f}');checks=[c for c in cc if c.get('reference')==r['id']];admitted=r['status']=='PASS' and len(checks)==2 and all(c['admitted'] for c in checks)
        freq.append((f'{f/1e6:g} MHz',r if admitted else None,'dynamic'))
    jobs=[('Video_F10AE_1_target_shape_optimization','T1: 3D aspect +3% +/-0.5 pp | COM <=20 nm', [('Initial state',t1,'initial'),('Actual accepted trajectory',t1,'dynamic'),('Actual endpoint',t1,'final')],.002),('Video_F10AE_2_deformation_and_recovery','Two forward actions | switch at 2 ms | quasi-steady radial response', [('Initial reference',t4,'initial'),('Active second action',t4,'dynamic'),('First action then OFF',by['schedule_C'],'dynamic')],.004),('Video_F10AE_3_frequency_and_robustness','Frozen T1 action | matched field and radial response at each frequency',freq,.002)]
    manifest=[]
    for name,title,panels,duration in jobs:
        begin=time.perf_counter();data=[traces(r) if r else None for label,r,mode in panels];times=np.arange(round(duration/5e-6)+1)*5e-6;path=out/(name+'.mp4');writer=av.open(str(path),'w');stream=writer.add_stream('libx264',rate=60);stream.width=1800;stream.height=1000;stream.pix_fmt='yuv420p';stream.options={'crf':'18','preset':'medium'};mapping=[]
        for frame,t in enumerate(times):
            im=Image.new('RGB',(1800,1000),'#080d16');d=ImageDraw.Draw(im);d.text((24,18),title,font=font(25),fill='#edf3fa');d.text((24,58),f'Physical time: {t*1e3:.3f} ms | N = 64 | R = 1.78 um | zero background flow | slow playback',font=font(21),fill='#c4d8eb');sources=[]
            for k,((label,r,mode),arr) in enumerate(zip(panels,data)):
                left=24+600*k;d.text((left,105),label,font=font(24),fill='#91c8f3')
                if arr is None:d.text((left,330),'NUMERICAL_NOT_ADMITTED',font=font(22),fill='#e0af7c');continue
                tt,x=arr;idx=0 if mode=='initial' else len(tt)-1 if mode=='final' else int(np.searchsorted(tt,t-1e-13))
                if mode=='dynamic':assert abs(tt[idx]-t)<1e-13
                m=measure(x[0],x[idx]);im.paste(scene(x[idx],size=(560,560)),(20+k*600,140));d=ImageDraw.Draw(im)
                stage=0
                segments=r['task'].get('segments')
                if segments:stage=int(tt[idx]>=segments[0]['duration']);action=segments[stage]
                else:action=r['task']
                gain=(m['aspect_3d']/r['initial']['aspect_3d']-1)*100;com=np.linalg.norm(m['center_3d_m']-r['initial']['center_3d_m'])*1e9;length=m['principal_rms_lengths_3d_m']*1e6
                lines=[f'Rg: {m["Rg_3d_m"]*1e6:.3f} um | aspect change: {gain:+.3f}%',f'L1/L2/L3: {length[0]:.3f} / {length[1]:.3f} / {length[2]:.3f} um',f'COM: {com:.2f} nm | Kabsch RMS: {m["kabsch_rmsd_3d_m"]*1e9:.2f} nm',f'Components: {m["proximity_components"]} | largest: {m["largest_component"]} | contact pairs: {m["hard_contact_pairs"]}', 'A [Y-/Y+/X-/X+]: '+' / '.join(f'{v:.0f}' for v in action['amplitudes'])+' Pa-eq','Phase: '+' / '.join(f'{np.degrees(v):.1f}' for v in action['phases'])+' deg',f'Stored sample: {idx} | t = {tt[idx]*1e3:.3f} ms']
                for j,line in enumerate(lines):d.text((left,708+j*29),line,font=font(17),fill='#d0dfeb')
                sources.append(dict(id=r['id'],trajectory=r['trajectory_path'],sha256=r['trajectory_sha256'],saved_index=idx,mode=mode))
            if name.endswith('optimization'):
                g=next(g for g in validation['groups'] if g['target']=='T1');status=f'Discovery: {"PASS" if t1["targets"]["T1_DISCOVERY_PASS"] else "FAIL"} | dt/cache: {"PASS" if g["replicates"][0]["numerical_pass"] else "FAIL"} | independent validation: {sum(r["confirmed"] for r in g["replicates"][1:])}/3 | confirmed control: {"PASS" if admission["targets"]["T1"]["confirmed"] else "FAIL"}'
            elif name.endswith('recovery'):
                status='Recovery criterion: '+('PASS at discovery; requires independent admission' if all(t4['recovery']['raw_recovery_pass'].values()) else 'NOT ACHIEVED in this trajectory')
            else:status='Only complete trajectories passing dt/cache comparisons are displayed.'
            d.text((24,929),status,font=font(20),fill='#e5bc89');d.text((24,969),'EXPLORATORY | Uncalibrated homogeneous water | Exact fixed-ID positions and sphere radii | training_ready=false',font=font(16),fill='#aebdcc')
            for packet in stream.encode(av.VideoFrame.from_image(im)):writer.mux(packet)
            if frame in [0,len(times)//2,len(times)-1]:im.save(out/f'{name}_{frame:04d}.png')
            mapping.append(dict(frame=frame,physical_time_s=float(t),sources=sources))
        for packet in stream.encode():writer.mux(packet)
        writer.close();pts=[]
        with av.open(str(path)) as dec:
            for fr in dec.decode(video=0):pts.append(float(fr.pts*fr.time_base))
        assert len(pts)==len(times) and np.allclose(pts,np.arange(len(times))/60,atol=1e-6)
        record=dict(name=name,sha256=sha(path),fps=60,frames=len(pts),physical_duration_s=duration,playback_s=len(pts)/60,resolution=[1800,1000],fully_decoded=True,uniform_PTS=True,uniform_physical_step_s=5e-6,mapping=mapping,wall_s=time.perf_counter()-begin,identity=identity())
        write(out/(name+'.json'),record);manifest.append({k:v for k,v in record.items() if k!='mapping'})
    write(out/'MANIFEST.json',manifest)
if __name__=='__main__':main()
