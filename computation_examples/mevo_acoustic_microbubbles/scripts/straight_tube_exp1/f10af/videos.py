"""Exact accepted HF positions and separately labeled LF dense output, 60 fps."""
import time
import numpy as np,h5py,av
from PIL import Image,ImageDraw,ImageFont
from .common import *
from .metrics import spatial,rotation
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
def font(n):return ImageFont.truetype(FONT,n)
def scene(x,size=(880,650)):
    # Orthographic XY projection. True physical radii, no cluster transformation.
    w,h=size;ss=2;scale=min(w/.00124,h/.00104);im=Image.new('RGB',(w*ss,h*ss),'#080d16');d=ImageDraw.Draw(im)
    def project(q):return ((w/2+scale*q[0])*ss,(h/2-scale*q[1])*ss)
    for xg in np.arange(-.0005,.00051,.0001):d.line([project([xg,-.0004]),project([xg,.0004])],fill='#25303f',width=1)
    for yg in np.arange(-.0004,.00041,.0001):d.line([project([-.0005,yg]),project([.0005,yg])],fill='#25303f',width=1)
    radius=R*scale*ss
    for q in x[np.argsort(x[:,2])]:
        px,py=project(q)
        # Draw actual-radius sphere footprint, with a small internal light highlight.
        d.ellipse((px-radius,py-radius,px+radius,py+radius),fill='#347dc0');rr=radius*.7;d.ellipse((px-rr-.15*radius,py-rr-.2*radius,px+rr-.15*radius,py+rr-.2*radius),fill='#63baff')
    c=x.mean(0);cx,cy=project(c);d.line((cx-5*ss,cy,cx+5*ss,cy),fill='#f2c784',width=ss);d.line((cx,cy-5*ss,cx,cy+5*ss),fill='#f2c784',width=ss)
    im=im.resize(size,Image.Resampling.LANCZOS);d=ImageDraw.Draw(im);d.line((30,h-28,30+scale*100e-6,h-28),fill='#e2eaf3',width=3);d.text((30,h-56),'100 µm',font=font(18),fill='#d2e0ed');d.text((w-320,h-48),'XY projection | true sphere radii',font=font(17),fill='#c1cedb');return im

def load(r,lf=False):
    if lf:
        with h5py.File(BASE/r['trajectory_path']) as h:return h['sample_time_s'][:],h['sample_positions_m'][:],None
    with h5py.File(BASE/'dataset'/(r['id']+'.h5')) as h:return h['time_s'][:],h['positions_m'][:],h['force_samples/velocity_m_s'][:]
def main():
    out=BASE/'videos';out.mkdir(exist_ok=True);rows=read(BASE/'ALL_HF_RUNS.json');good=[r for r in rows if r['status']=='PASS' and r['task']['state']=='N64_wide_r0' and not r['task'].get('variant') and not r['task'].get('offset_m')];by={r['id']:r for r in rows}
    def matching(a,b):
        choices=[(r,s) for r in good for s in good if r['task'].get('action')==a and s['task'].get('action')==b and r['completed_s']==s['completed_s']];return max(choices,key=lambda t:t[0]['completed_s'])
    p1=matching('X_PAIR','XY_QUADRATURE');p2=matching('ROT_PLUS','ROT_MINUS');x=max([r for r in good if r['task'].get('action')=='X_PAIR' and (BASE/'lf'/r['id']/'SUMMARY.json').exists() and read(BASE/'lf'/r['id']/'SUMMARY.json')['status']=='PASS'],key=lambda r:r['completed_s']);lf=read(BASE/'lf'/x['id']/'SUMMARY.json')
    sci=read(BASE/'SCIENTIFIC_RESULTS.json');science_by={r['id']:r for r in sci['rows']}
    jobs=[('Video_F10AF_1_stripe_and_pattern_formation','Actual pattern-formation attempts',[(p1[0],False,'X pair'),(p1[1],False,'XY quadrature')]),('Video_F10AF_2_collective_rotation_feasibility','Opposite spatial phase orders',[(p2[0],False,'Positive phase order'),(p2[1],False,'Negative phase order')]),('Video_F10AF_3_multifidelity_vs_full_physics','Same action, initial state and physical time',[(x,False,'Full F10AA physics'),(lf,True,'LOW_FIDELITY_PREDICTION')])];manifest=[]
    for name,title,panels in jobs:
        begin=time.perf_counter();data=[load(r,l) for r,l,label in panels];T=min(q[0][-1] for q in data);tt=np.linspace(0,T,401);path=out/(name+'.mp4');writer=av.open(str(path),'w');stream=writer.add_stream('libx264',rate=60);stream.width=1920;stream.height=1080;stream.pix_fmt='yuv420p';stream.options={'crf':'18','preset':'medium'};mapping=[];initial=[spatial(d[1][0]) for d in data]
        for frame,t in enumerate(tt):
            im=Image.new('RGB',(1920,1080),'#080d16');d=ImageDraw.Draw(im);d.text((35,20),title,font=font(31),fill='#edf3fa');d.text((35,65),f'Physical time: {t*1e3:.3f} ms   |   2.5 MHz   |   N = 64   |   R = 1.78 µm   |   zero background flow',font=font(23),fill='#c5d8eb');sources=[];positions=[];panel_metrics=[]
            for k,((r,islf,label),(times,xx,vv)) in enumerate(zip(panels,data)):
                j=int(np.argmin(abs(times-t)));assert abs(times[j]-t)<1e-12,'Uniform movie times must correspond to saved positions';xj=xx[j];positions.append(xj);m=spatial(xj);panel_metrics.append(m);left=35+950*k;im.paste(scene(xj),(left,150));d=ImageDraw.Draw(im);d.text((left,115),label,font=font(27),fill='#8dc9f5');task=r['task'];a=' / '.join(f'{v:g}' for v in task['amplitudes']);p=' / '.join(f'{np.degrees(v):g}' for v in task['phases']);com=np.linalg.norm(xj.mean(0)-xx[0].mean(0))*1e6
                lines=[f'Y- / Y+ / X- / X+ amplitudes: {a} Pa-eq',f'Phases: {p} deg',f'COM shift: {com:.3f} um | Rg: {m["Rg_3d_m"]*1e6:.2f} um',f'Peak S: {m["structure_peak_S"]:.3f} | candidate spacing: {m["stripe_spacing_m"]*1e6:.1f} um']
                if name.endswith('rotation_feasibility'):
                    rot=rotation(xj,vv[j]);lines.append(f'Omega_z: {rot["Omega_z_s"]:+.4f} rad/s | rigid fraction: {rot["rigid_fraction"]:.4f}')
                elif name.endswith('full_physics'):lines[3]=f'Minimum surface gap: {m["minimum_gap_m"]*1e9:.3f} nm | numerical status: PASS';lines.append('Stored HF accepted states' if not islf else 'LF dense output; original solver near contact')
                else:lines.append('Numerical status: PASS | Stored HF accepted states' if not islf else 'SciPy RK45 dense output / labeled HF-fallback interpolation')
                for n,line in enumerate(lines):d.text((left,815+32*n),line,font=font(19),fill='#d4e1ed')
                sources.append(dict(id=r['id'],LF=islf,source_index=j,source_time_s=float(times[j]),raw_sha256=r['trajectory_sha256']))
            if name.endswith('full_physics'):status=f'Maximum fixed-ID difference: {np.linalg.norm(positions[0]-positions[1],axis=1).max()*1e9:.3f} nm | Rg difference: {abs(panel_metrics[0]['Rg_3d_m']-panel_metrics[1]['Rg_3d_m'])*1e9:.3f} nm | HF {x["task_wall_s"]:.2f} s; LF {lf["wall_s"]:.2f} s'
            elif name.endswith('rotation_feasibility'):status=('Rotation candidate requires independent confirmation' if any(science_by[r['id']]['rotation_candidate'] for r,_,_ in panels) else 'COLLECTIVE_ROTATION_NOT_DEMONSTRATED')+' | COM markers are reference centers.'
            else:status=('Pattern candidate requires independent confirmation.' if any(science_by[r['id']]['stripe_candidate'] or science_by[r['id']]['lattice_candidate'] for r,_,_ in panels) else 'Stripe / lattice formation not confirmed by registered criteria.')+' Actual simulated motion.'
            d.text((35,1000),status,font=font(21),fill='#e3bd8c');d.text((35,1045),'EXPLORATORY | Homogeneous-water proxy | No streaming or spin | Real camera observability unverified | training_ready=false',font=font(17),fill='#b7c5d3')
            for packet in stream.encode(av.VideoFrame.from_image(im)):writer.mux(packet)
            if frame in [0,200,400]:im.save(out/f'{name}_{frame:04d}.png')
            mapping.append(dict(frame=frame,physical_time_s=float(t),sources=sources))
        for packet in stream.encode():writer.mux(packet)
        writer.close();pts=[]
        with av.open(str(path)) as reader:
            for f in reader.decode(video=0):pts.append(float(f.pts*f.time_base))
        assert len(pts)==401 and np.allclose(pts,np.arange(401)/60,atol=1e-7)
        r=dict(file=path.name,sha256=sha(path),frames=401,fps=60,resolution=[1920,1080],physical_duration_s=float(T),uniform_physical_step_s=float(T/400),playback_duration_s=401/60,fully_decoded=True,uniform_PTS=True,rendering='Exact particle coordinates, true radius 1.78 µm, fixed XY projection; no clustering simplification or coordinate rotation',mapping=mapping,wall_s=time.perf_counter()-begin,identity=identity());write(out/(name+'.json'),r);manifest.append({k:v for k,v in r.items() if k!='mapping'});print('VIDEO',name,flush=True)
    write(out/'MANIFEST.json',dict(videos=manifest,identity=identity()))
if __name__=='__main__':main()
