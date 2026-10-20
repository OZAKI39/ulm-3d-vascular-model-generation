"""Accepted native Y-clot trajectories only; red/Z-up publication rendering."""
from common import *
from current_surface import CurrentSurface
import sys,time,csv,subprocess,argparse,html
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageDraw

# Reuse the last approved red presentation, not its old paths or its states.
red=module('approved_red_visual_style',R/'analysis/clot_fragment_surface_coupling_v3_red_visuals/render.py')
graph=module('y_application_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
CASE=O/'cases/biphasic';VIS=F/'figures';VIDEO=F/'videos'

def catalog():
    found={}
    for path in sorted((CASE/'keys').glob('*.npz')):
        if path.stem=='failed_unaccepted' or '.tmp.' in path.name:continue
        with np.load(path) as z:step=int(z['step']);t=float(z['time'])
        if step not in found or path.stem in ['initial','final','last_accepted','peak_load']:
            found[step]=dict(path=path,time=t,step=step)
    return sorted(found.values(),key=lambda q:q['time'])

def load(path):
    with np.load(path) as z:return {k:z[k].copy() for k in red.old.FIELDS+['surface_traction','wall_clearance'] if k in z}

def save_image(name,array):
    im=Image.fromarray(array);im.save(VIS/(name+'.png'));im.convert('RGB').save(VIS/(name+'.pdf'),resolution=240.)

def draw(ref,surf,item,bounds,field=False,size=(1920,1080)):
    s=load(item['path']);geom=surf.geometry(s['x'],s['d'])
    assert np.array_equal(geom['active'],s['surface_active'])
    assert np.allclose(geom['area'][geom['active']],s['surface_area'][geom['active']],rtol=1e-12,atol=1e-24)
    pl=red.setup(bounds,size)
    wall=surf.wall.wall.copy();wall.points*=1e6
    clipped=wall.clip_box(bounds,invert=False).extract_surface()
    if clipped.n_points:
        pl.add_mesh(clipped,color='#aabdc7',opacity=.12,lighting=False)
    stats=red.scene(pl,ref,surf,s,geom,field=field,arrow_limit=60 if field else 38)
    pl.render();raw=pl.screenshot(return_img=True);picture=red.annotate(pl,raw,bounds,s,stats,field)
    stats.update(source=str(item['path'].relative_to(R)),source_sha256=sha(item['path']),physical_coordinates_unscaled=True)
    pl.close();return picture,stats

def overview(ref,surf):
    wall=surf.wall.wall.copy();wall.points*=1e6
    bounds=np.array(wall.bounds);center=np.array(wall.center);span=max(np.ptp(wall.points,axis=0))
    pl=red.setup(bounds,size=(2400,1500));pl.add_mesh(wall,color='#b9cbd5',opacity=.17,smooth_shading=True)
    points=ref['X'][:int(ref['clot_count'])]*1e6
    pl.add_points(points,color=red.RED,point_size=7,render_points_as_spheres=True)
    direction=np.array([.8,-1.6,2.4]);pl.camera_position=[center+span*direction,center,(0,0,1)]
    forward=-direction/np.linalg.norm(direction);right=np.cross(forward,[0.,0.,1.]);right/=np.linalg.norm(right);up=np.cross(right,forward)
    corners=np.array([(x,y,z) for x in bounds[:2] for y in bounds[2:4] for z in bounds[4:]])
    aspect=(.905*2400)/(.78*1500);pl.camera.parallel_scale=1.16*max(np.ptp(corners@up)/2,np.ptp(corners@right)/(2*aspect))
    pl.render();raw=pl.screenshot(return_img=True)
    im=Image.new('RGB',(2400,1500),'white');view=(132,135,2304,1305);im.paste(Image.fromarray(raw).crop(view),view[:2]);ink=ImageDraw.Draw(im)
    ink.text((60,40),'Original Y vessel and irregular wall-attached clot',font=red.font(36),fill='#28333b')
    point=red.project(pl,points.mean(0));label=point+[-200,-115]
    ink.line([tuple(label+[50,35]),tuple(point)],fill=red.RED,width=3);ink.text(tuple(label),'Clot',font=red.font(34),fill=red.RED)
    lo,hi=bounds[::2],bounds[1::2]
    for axis in range(3):
        first=lo.copy();last=lo.copy()
        if axis==1:first[0]=last[0]=hi[0]
        last[axis]=hi[axis]
        anchor=['mt','lm','rm'][axis];shift=np.array([(0,18),(23,7),(-23,0)][axis])
        ticks=np.linspace(lo[axis],hi[axis],3 if axis==2 else 4)
        for k,value in enumerate(ticks):
            if axis in (1,2) and k==0:continue
            q=first.copy();q[axis]=value;pos=red.project(pl,q)+shift
            ink.text(tuple(pos),f'{value/1000:.1f}',font=red.font(30),anchor=anchor,fill='#34454f')
        extra=np.array([(0,60),(90,50),(-100,0)][axis]);pos=red.project(pl,(first+last)/2)+extra
        ink.text(tuple(pos),'XYZ'[axis]+' (mm)',font=red.font(32),anchor=anchor,fill='#24333c')
    p0=red.project(pl,center);p1=red.project(pl,center+[0,0,1000]);assert abs(p0[0]-p1[0])<1e-7
    pl.close();save_image('FIG_01_Y_VESSEL_INITIAL_GEOMETRY',np.asarray(im))

def curves():
    with (CASE/'history.csv').open() as f:rows=list(csv.DictReader(f))
    def a(k):return np.array([float(q[k]) for q in rows])
    t=a('time_s')*1e3
    fig,axs=plt.subplots(2,3,figsize=(16,8),constrained_layout=True)
    ax=axs[0,0];ax.plot(t,a('nominal_shear_Pa'),label='Prescribed nominal');ax.plot(t,a('shear_peak_Pa'),label='Actual surface peak');ax.set_ylabel('Shear (Pa)');ax.legend()
    ax=axs[0,1];ax.plot(t,a('bulk_failed'),label='Bulk');ax.plot(t,a('interface_failed'),label='Attachment');ax.set_ylabel('Failed original bonds');ax.legend()
    ax=axs[0,2];ax.plot(t,a('detached_nodes')/1927*100,color=red.RED);ax.set_ylabel('Detached material (%)')
    ax=axs[1,0];ax.plot(t,a('largest_nodes'),label='Largest free component');ax.set_ylabel('Original nodes');ax.legend()
    ax=axs[1,1];ax.plot(t,a('N20_200'),color=red.RED);ax.set_ylabel('Free components: 20–200 nodes')
    ax=axs[1,2];ax.plot(t,a('relative_R')*100,color='#17669a');ax.set_ylabel('Energy residual (%)',color='#17669a')
    other=ax.twinx();other.plot(t,a('L2_active_patches'),color=red.RED);other.set_ylabel('Active L2 patches',color=red.RED)
    for ax in axs.ravel():ax.set_xlabel('Physical time (ms)');ax.grid(alpha=.3)
    fig.savefig(VIS/'FIG_05_DYNAMICS_AND_AUDIT.png',dpi=200);fig.savefig(VIS/'FIG_05_DYNAMICS_AND_AUDIT.pdf');plt.close(fig)

def main(preview=False):
    start=time.perf_counter();VIS.mkdir(parents=True,exist_ok=True);VIDEO.mkdir(exist_ok=True)
    ref=load_npz(O/'provenance/reference.npz');surf=CurrentSurface(ref,'uniform');items=catalog();assert items
    allpoints=[ref['X'][int(ref['clot_count']):]]
    for q in items:
        with np.load(q['path']) as z:allpoints.append(z['x'][:int(ref['clot_count'])])
    bounds=red.bounds_for(np.vstack(allpoints)*1e6,pad=.10)
    if preview:
        pic,stats=draw(ref,surf,items[0],bounds)
        Image.fromarray(pic).save(VIS/'PREVIEW_INITIAL.png');overview(ref,surf)
        dump(O/'PREVIEW_RECORD.json',stats);return
    summary=json.loads((CASE/'SUMMARY.json').read_text())
    peak=min(items,key=lambda q:abs(q['time']-HOLD));mid=min(items,key=lambda q:abs(q['time']-OFF))
    records=[];overview(ref,surf)
    for name,item,field in [('FIG_02_INITIAL_CLOT',items[0],False),('FIG_03_PEAK_LOAD_FIELD',peak,True),('FIG_04_FINAL_CLOT',items[-1],False)]:
        pic,stats=draw(ref,surf,item,bounds,field=field,size=(3200,1800));save_image(name,pic);stats['figure']=name;records.append(stats)
        if field:
            s=load(item['path']);ids=np.flatnonzero(s['surface_active'])
            with (F/'PEAK_SURFACE_LOADS.csv').open('w') as f:
                w=csv.writer(f);w.writerow(['face_id','method','x_m','y_m','z_m','area_m2','opening_weight','tx_Pa','ty_Pa','tz_Pa','Fx_N','Fy_N','Fz_N'])
                for i in ids:w.writerow([i,int(surf.method[i]),*s['surface_centroid'][i],s['surface_area'][i],s['surface_weight'][i],*s['surface_traction'][i],*s['surface_force'][i]])
            err=float(np.max(abs(s['surface_traction']*s['surface_area'][:,None]-s['surface_force'])))
            assert err<1e-16;stats['traction_times_geometric_area_error_N']=err
    curves()
    # Real saved frames, one per state, and an explicit final hold. No invented
    # intermediate material state, no time interpolation, no position scaling.
    order=items+[items[-1]]*36;target=VIDEO/'VIDEO_Y_CLOT_CURRENT_METHOD_RED_Z_VERTICAL.mp4'
    cmd=['ffmpeg','-y','-loglevel','error','-f','rawvideo','-vcodec','rawvideo','-pix_fmt','rgb24','-s','1920x1080','-r','24','-i','-',
         '-an','-c:v','libx264','-crf','18','-preset','medium','-pix_fmt','yuv420p','-movflags','+faststart',str(target)]
    pipe=subprocess.Popen(cmd,stdin=subprocess.PIPE);video_records=[];last=None;cached=None
    try:
        for i,item in enumerate(order):
            if item['path']!=last:cached,stats=draw(ref,surf,item,bounds);last=item['path']
            pipe.stdin.write(cached.tobytes());video_records.append(dict(frame=i,time_s=item['time'],source=str(item['path'].relative_to(R)),source_sha256=sha(item['path']),z_axis_horizontal_error_px=stats['z_axis_horizontal_error_px']))
            if i in [0,len(items)//2,len(items)-1]:Image.fromarray(cached).save(VIS/f'VIDEO_FRAME_{i:03d}.png')
        pipe.stdin.close();assert pipe.wait()==0
    finally:
        if pipe.poll() is None:pipe.kill()
    subprocess.run(['ffmpeg','-v','error','-i',str(target),'-f','null','-'],check=True,capture_output=True)
    dump(O/'VISUALIZATION_VALIDATION.json',dict(status='PASS',formal_status=summary['status'],accepted_end_time_s=summary['time_s'],
        frames=len(order),unique_saved_states=len(items),resolution=[1920,1080],fps=24,frame_sources=video_records,static_figures=records,
        all_material_nodes_visible=True,drawn_bonds='Surviving original short bonds only; full graph used for topology',Z_vertical=True,
        displacement_scale=1.,interpolated_physics_frames=0,pdf='Figures 1–4 embed raster; dynamics plot is vector',wall_s=time.perf_counter()-start))
    links=''.join(f'<h2>{html.escape(p.stem)}</h2><a href="figures/{p.name}"><img src="figures/{p.name}" style="max-width:100%"></a>' for p in sorted(VIS.glob('FIG_*.png')))
    text=f'''<!doctype html><html lang="zh"><meta charset="utf-8"><title>Y 血管原血块：当前方法应用</title>
<body style="max-width:1300px;margin:auto;font-family:sans-serif;line-height:1.7"><h1>Y 形血管内原不规则血块：当前方法应用</h1>
<p>状态：{summary['status']}；接受物理时间：{summary['time_s']*1e3:.6f} ms。红色点云和存活材料键均来自原生接受态；Z 轴竖直、位移 ×1。</p>
<p>载荷为规定的 1000 Pa 机械剪切代理。原 Y 管腔与壁接触参与；没有重新求解 CFD、微泡声流或双向 FSI。强区为原场平移裁剪，比例 45.7706%。</p>
<video controls style="width:100%" src="videos/{target.name}"></video>{links}
<p><a href="PEAK_SURFACE_LOADS.csv">实际峰值阶段表面载荷 CSV</a></p></body></html>'''
    (F/'OPEN_RESULTS.html').write_text(text)
    print('Y_VISUALS_COMPLETE',target,len(order),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preview',action='store_true');args=p.parse_args();main(args.preview)
