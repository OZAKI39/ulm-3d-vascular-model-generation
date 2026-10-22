"""Actual fixed-ID trajectories; orthographic true-radius sphere rendering, no interpolation."""
import time,math
import numpy as np,h5py,av
from PIL import Image,ImageDraw,ImageFont
from .common import *
from mevo_mb.fast_sim.f10g_contact import resistance
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'

def font(size):return ImageFont.truetype(FONT,size)
def load_run(backend,n=64):
    p=BASE/f'full/{backend}/runs/N{n}_PAIR_2000us_DT1'
    with h5py.File(p/'trajectory.h5') as h:
        return dict(path=p,x=h['positions_m'][:],t=h['time_s'][:],metrics=h['step_metrics'][:],columns=h['step_metrics'].attrs['columns'].split(','),force=h['force_samples/total_force_N'][:],contact=h['force_samples/pair_contact_reaction_N'][:])

def scene(x,size=(590,590),contact=None):
    width,height=size;im=np.zeros((height,width,3),dtype=np.uint8)+np.array([8,13,22],dtype=np.uint8);depth=np.full((height,width),-np.inf)
    az,el=np.deg2rad([-60,27]);view=np.array([np.cos(el)*np.cos(az),np.cos(el)*np.sin(az),np.sin(el)]);right=np.array([-np.sin(az),np.cos(az),0]);up=np.cross(view,right)
    transform=np.stack([right,up,view]);scale=min(size)/105e-6;q=x@transform.T;centers=np.column_stack([width/2+scale*q[:,0],height/2-scale*q[:,1]])
    image=Image.fromarray(im);d=ImageDraw.Draw(image)
    # Fixed 20 um grid, identical camera, viewport and physical sphere radius in both panels.
    for axis in [0,1]:
        for value in [-40e-6,-20e-6,0,20e-6,40e-6]:
            pts=[]
            for end in [-40e-6,40e-6]:
                p=np.array([value,end,-32e-6]) if axis==0 else np.array([end,value,-32e-6]);pp=transform@p;pts.append((width/2+scale*pp[0],height/2-scale*pp[1]))
            d.line(pts,fill=(35,49,66),width=1)
    im=np.array(image)
    if contact is not None:
        d=ImageDraw.Draw(image);ii,jj=np.triu_indices(len(x),1)
        for i,j in zip(ii[contact>0],jj[contact>0]):d.line([tuple(centers[i]),tuple(centers[j])],fill=(238,187,87),width=2)
        im=np.array(image)
    radius=R*scale;light=np.array([-.35,.5,.79]);light/=np.linalg.norm(light)
    for k,(cx,cy) in enumerate(centers):
        left=max(0,int(np.floor(cx-radius-1)));rightp=min(width,int(np.ceil(cx+radius+1)));top=max(0,int(np.floor(cy-radius-1)));bottom=min(height,int(np.ceil(cy+radius+1)))
        yy,xx=np.mgrid[top:bottom,left:rightp];nx=(xx+.5-cx)/radius;ny=-(yy+.5-cy)/radius;r2=nx*nx+ny*ny;inside=r2<=1;nz=np.sqrt(np.maximum(0,1-r2));z=q[k,2]+R*nz
        visible=inside&(z>depth[top:bottom,left:rightp]);diff=np.maximum(0,nx*light[0]+ny*light[1]+nz*light[2]);spec=np.maximum(0,.3*nx+.45*ny+.84*nz)**20
        rgb=np.clip(np.array([42,132,219])*(.42+.58*diff[:,:,None])+95*spec[:,:,None],0,255).astype(np.uint8)
        im[top:bottom,left:rightp][visible]=rgb[visible];depth[top:bottom,left:rightp][visible]=z[visible]
    image=Image.fromarray(im);d=ImageDraw.Draw(image)
    if contact is not None:
        ii,jj=np.triu_indices(len(x),1)
        for i,j in zip(ii[contact>0],jj[contact>0]):d.line([tuple(centers[i]),tuple(centers[j])],fill=(238,187,87),width=2)
    d.line((30,height-36,30+20e-6*scale,height-36),fill='white',width=3);d.text((30,height-65),'20 um',font=font(17),fill='#dce5ee');d.text((width-190,height-45),'3D orthographic',font=font(16),fill='#b8c6d4')
    return image

def chart(d,box,t,values,current,title,unit,color):
    x0,y0,x1,y1=box;d.text((x0,y0-27),title,font=font(18),fill='#e0e7ef');d.line((x0,y0,x0,y1,x1,y1),fill='#76889b',width=1)
    lo=min(0,float(np.min(values)));hi=float(np.max(values));hi=hi if hi>lo else lo+1
    pts=[(x0+(x1-x0)*a/t[-1],y1-(y1-y0)*(v-lo)/(hi-lo)) for a,v in zip(t,values)];d.line(pts,fill=color,width=2)
    xx=x0+(x1-x0)*current/t[-1];d.line((xx,y0,xx,y1),fill='#cccccc',width=1)
    d.text((x0,y1+4),'0',font=font(14),fill='#aebdcb');d.text((x1,y1+4),f'{t[-1]*1e3:.1f} ms',font=font(14),fill='#aebdcb',anchor='ra');d.text((x1-100,y0),f'{hi:.2g} {unit}',font=font(14),fill='#c9d5e1')

def main():
    out=BASE/'videos';out.mkdir(exist_ok=True);ref=load_run('B0');opt=load_run('B1');n=64;indices=np.arange(0,len(opt['t']),10);fps=30;records=[]
    for mode,name in [(1,'Video_F10AA_1_reference_vs_accelerated'),(2,'Video_F10AA_2_contact_and_solver_diagnostics')]:
        p=out/(name+'.mp4');begin=time.perf_counter();writer=av.open(str(p),'w');stream=writer.add_stream('libx264',rate=fps);stream.width=1280;stream.height=800;stream.pix_fmt='yuv420p';stream.options={'crf':'18','preset':'medium'}
        for index in indices:
            im=Image.new('RGB',(1280,800),'#080d16');d=ImageDraw.Draw(im);t=opt['t'][index];k=min(max(0,index-1),len(opt['metrics'])-1)
            d.text((30,20),'F10AA | 64 bubbles | 2.5 MHz | zero background flow',font=font(24),fill='#ecf3fa');d.text((30,56),f'Physical time: {t*1e3:.3f} ms   |   N = 64   |   R = 1.78 um',font=font(20),fill='#c8d9e8')
            if mode==1:
                im.paste(scene(ref['x'][index]),(25,125));im.paste(scene(opt['x'][index]),(665,125));d=ImageDraw.Draw(im)
                d.text((30,95),'Reference CPU',font=font(22),fill='#d6dde5');d.text((670,95),'Accelerated matrix-free CPU',font=font(22),fill='#70b8ef')
                err=np.linalg.norm(ref['x'][index]-opt['x'][index],axis=1).max()*1e9
                d.text((30,727),f'Fixed-ID maximum difference: {err:.3g} nm | Same camera and physical radii',font=font(18),fill='#d6e2ec')
            else:
                im.paste(scene(opt['x'][index],contact=opt['contact'][k]),(20,125));d=ImageDraw.Draw(im)
                col=opt['columns'];metrics=opt['metrics'];tt=opt['t'][1:]
                chart(d,(700,150,1220,240),tt,metrics[:,col.index('active_count')],t,'Active hard-contact constraints','', '#e6bc65')
                chart(d,(700,325,1220,415),tt,metrics[:,col.index('force_residual_N')],t,'Force-balance residual','N','#68b4e9')
                chart(d,(700,500,1220,590),tt,metrics[:,col.index('actual_complementarity_m')]*1e9,t,'True-gap complementarity','nm','#77c1b0')
                force=np.linalg.norm(opt['force'][k],axis=1).max();active=int(metrics[k,col.index('active_count')]);refms=ref['metrics'][k,ref['columns'].index('full_step_s')]*1e3;fastms=metrics[k,col.index('full_step_s')]*1e3
                d.text((30,97),'Orange links: active pair constraints',font=font(18),fill='#e6bc65')
                drag=np.trace(resistance(opt['x'][k],R)['A'])/(3*n)
                d.text((700,718),f'Resistance: Tr(A)/(3N) = {drag:.3g}',font=font(17),fill='#d6e2ec')
                d.text((30,725),f'Active pairs: {active} | Maximum acoustic force: {force*1e12:.3g} pN',font=font(18),fill='#d6e2ec')
                d.text((700,660),f'Observed step: CPU reference {refms:.2f} ms',font=font(17),fill='#d6e2ec');d.text((700,688),f'Optimized CPU {fastms:.2f} ms',font=font(17),fill='#70b8ef')
            d.text((30,766),'EXPLORATORY | Uncalibrated homogeneous-water proxy | 10 us per displayed frame | training_ready=false',font=font(15),fill='#a8b6c5')
            frame=av.VideoFrame.from_ndarray(np.array(im),format='rgb24')
            for packet in stream.encode(frame):writer.mux(packet)
            if index in [indices[0],indices[len(indices)//2],indices[-1]]:im.save(out/f'{name}_{index:04d}.png')
        for packet in stream.encode():writer.mux(packet)
        writer.close();count=0;pts=[]
        with av.open(str(p)) as reader:
            for frame in reader.decode(video=0):count+=1;pts.append(float(frame.pts*frame.time_base))
        assert count==len(indices) and np.allclose(np.diff(pts),1/fps,atol=1e-9)
        records.append(dict(file=p.name,frames=count,fps=fps,playback_duration_s=count/fps,physical_start_s=0,physical_end_s=float(opt['t'][-1]),uniform_physical_step_s=1e-5,source_frames=indices.tolist(),sha256=sha(p),decode_all_frames='PASS',strictly_uniform_PTS=True,wall_s=time.perf_counter()-begin,rendering='Actual per-ID positions and true 1.78um spheres; fixed orthographic camera and lighting, no interpolation or cluster rescaling'))
    write(out/'MANIFEST.json',dict(videos=records,reference_sha256=sha(ref['path']/'trajectory.h5'),accelerated_sha256=sha(opt['path']/'trajectory.h5'),source_sha256=sha(Path(__file__)),identity=identity()))
if __name__=='__main__':main()
