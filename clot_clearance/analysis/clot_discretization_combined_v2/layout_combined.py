"""Recompose the existing scientific scenes with larger editable labels.

Extracts the untouched, unannotated 3D scene embedded in each existing PDF.
No geometry or camera changes, no new scientific or 3D rendering run.
"""
from pathlib import Path
from pd_clot.portable_paths import relocate
import argparse,hashlib,io,json
import numpy as np
from PIL import Image
import pymupdf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'results/clot_discretization_three_views_v1'
OUT=ROOT/'results/clot_discretization_combined_v2'
FONT=Path(str(relocate('/mnt/c/Windows/Fonts/arial.ttf')))
FIGURE_MM=(540,195)
SIZES=dict(title=28,view_title=24,lattice_label=22,ticks=18,axis_labels=20,legend=20)
COLORS=dict(clot='#b85c64',anchor='#8aa7b9',debris='#586d7b',bond='#718593',text='#263746')

def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def scene_rgb(doc,xref):
    data=doc.extract_image(xref);arr=np.asarray(Image.open(io.BytesIO(data['image'])).convert('RGB'))
    info=next(r for r in doc[0].get_image_info(xrefs=True) if r['xref']==xref);a,b,c,d,_,_=info['transform']
    assert abs(b)<1e-10 and abs(c)<1e-10
    # Matplotlib stores its PDF scene upside down and uses a negative-height
    # image transform. Restore the displayed orientation before recomposition.
    if d<0:arr=arr[::-1]
    if a<0:arr=arr[:,::-1]
    return np.ascontiguousarray(arr)

def project(point,a,box):
    delta=np.asarray(point)-a['camera_focal_mm'];aspect=a['size'][0]/a['size'][1]
    u=.5+delta@a['camera_right']/(2*a['parallel_scale_mm']*aspect)
    v=.5+delta@a['camera_up']/(2*a['parallel_scale_mm'])
    x,y,w,h=box;return np.array([x+w*u,y+h*v])

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--preview',action='store_true');args=parser.parse_args();dpi=120 if args.preview else 600
    OUT.mkdir(parents=True,exist_ok=True);audits=json.loads((SOURCE/'PROVENANCE.json').read_text())['views'];sources=[SOURCE/'PROVENANCE.json',*[SOURCE/f'Fig_clot_discretization_view_{i}.pdf' for i in [1,2,3]]];before={str(p):sha(p) for p in sources}
    if FONT.exists():font_manager.fontManager.addfont(str(FONT))
    plt.rcParams.update({'font.family':'Arial' if FONT.exists() else 'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none','font.size':18,'text.color':'#20262d'})
    fig=plt.figure(figsize=np.asarray(FIGURE_MM)/25.4,dpi=dpi,facecolor='white');texts=[];pixel_hashes=[];panels=[]
    texts.append(fig.text(.016,.978,'Clot discretization',fontsize=SIZES['title'],va='top'))
    texts.append(fig.text(.984,.978,'Δx = 50 µm',fontsize=SIZES['lattice_label'],va='top',ha='right',color='#425665'))
    ticks=[[-6.8,-6.0,-5.2],[-.5,0,.5],[.6,1.,1.4]]
    for j,a in enumerate(audits):
        doc=pymupdf.open(sources[j+1]);images=doc[0].get_images(full=True);assert len(images)==1
        arr=scene_rgb(doc,images[0][0]);pixel_hashes.append(hashlib.sha256(arr.tobytes()).hexdigest())
        physical_width_mm=180*.94;height_mm=physical_width_mm*arr.shape[0]/arr.shape[1]
        box=[(j+.025)/3,.145,.94/3,height_mm/FIGURE_MM[1]]
        ax=fig.add_axes(box);ax.imshow(arr,interpolation='antialiased' if args.preview else 'none',aspect='equal');ax.axis('off')
        texts.append(fig.text((j+.5)/3,.895,f'View {j+1}',fontsize=SIZES['view_title'],ha='center',va='center',color=COLORS['text']))
        for item in a['axis_edges']:
            k=item['axis'];normal=np.asarray(item['outward_screen']);base=np.mean(item['edge'],axis=0)
            for value in ticks[k]:
                point=base.copy();point[k]=value;label=f'{value:.1f}' if k==0 else f'{value:g}'
                texts.append(ax.annotate(label,xy=project(point,a,box),xycoords=fig.transFigure,xytext=normal*22,textcoords='offset points',fontsize=SIZES['ticks'],color=COLORS['text'],ha='center',va='center',annotation_clip=False))
            if k==2:
                # Place the larger Z title above its axis rather than letting it
                # collide with the next panel or the outer page edge.
                endpoint=max(item['edge'],key=lambda q:project(q,a,box)[1]);anchor=project(endpoint,a,box);offset=(8,25)
            else:anchor=project(base,a,box);offset=normal*65
            texts.append(ax.annotate('XYZ'[k]+' (mm)',xy=anchor,xycoords=fig.transFigure,xytext=offset,textcoords='offset points',fontsize=SIZES['axis_labels'],color=COLORS['text'],ha='center',va='center',annotation_clip=False))
        panels.append(dict(view=j+1,box=box,source_pdf=str(sources[j+1]),scene_shape=list(arr.shape),scene_pixel_sha256=pixel_hashes[-1],physical_scene_width_mm=physical_width_mm,physical_scene_height_mm=height_mm))
    handles=[Line2D([],[],marker='o',linestyle='none',markersize=8,color=COLORS['clot'],label='Clot nodes'),Line2D([],[],marker='o',linestyle='none',markersize=8,color=COLORS['anchor'],label='Fixed anchors'),Line2D([],[],marker='o',linestyle='none',markersize=9,color=COLORS['debris'],label='Detached nodes'),Line2D([],[],color=COLORS['bond'],lw=1.6,label='Surviving bonds')]
    legend=fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.025),ncol=4,frameon=False,fontsize=SIZES['legend'],handlelength=1.1,columnspacing=1.8,handletextpad=.55)
    assert len(fig.legends)==1;texts.extend(legend.get_texts());fig.canvas.draw();renderer=fig.canvas.get_renderer();width,height=fig.canvas.get_width_height()
    boxes=[]
    for art in texts:
        b=art.get_window_extent(renderer);assert b.x0>=0 and b.y0>=0 and b.x1<=width+1 and b.y1<=height+1,(art.get_text(),b);boxes.append((art.get_text(),b))
    overlaps=[]
    for i,(t,b) in enumerate(boxes):
        for u,c in boxes[i+1:]:
            if min(b.x1,c.x1)-max(b.x0,c.x0)>2 and min(b.y1,c.y1)-max(b.y0,c.y0)>2:overlaps.append([t,u])
    assert not overlaps,overlaps
    stem='Fig_clot_discretization_three_views';name=stem+('_preview' if args.preview else '')
    fig.savefig(OUT/(name+'.png'),dpi=dpi)
    if not args.preview:
        fig.savefig(OUT/(stem+'.pdf'),dpi=dpi,metadata={'Title':'Clot discretization: three views with one shared legend'});fig.savefig(OUT/(stem+'.svg'),dpi=dpi)
        im=Image.open(OUT/(stem+'.png')).convert('RGB');im.thumbnail((2880,1600));im.save(OUT/(stem+'_preview.png'))
        # PDF export embeds the three exact original RGB scene arrays.
        doc=pymupdf.open(OUT/(stem+'.pdf'));out_hashes=[]
        for image_info in doc[0].get_images(full=True):
            rgb=scene_rgb(doc,image_info[0]);out_hashes.append(hashlib.sha256(rgb.tobytes()).hexdigest())
        assert sorted(out_hashes)==sorted(pixel_hashes)
        text=doc[0].get_text();assert all(text.count(label)==1 for label in ['Clot nodes','Fixed anchors','Detached nodes','Surviving bonds'])
    plt.close(fig);assert before=={p:sha(p) for p in before}
    audit=dict(status='PASS',source_sha256=before,figure_mm=FIGURE_MM,dpi=dpi,font_sizes_pt=SIZES,previous_font_sizes_pt=dict(title=12,ticks=8,axis_labels=9,legend=8),shared_legend_count=1,legend_position='Centered beneath middle panel; left and right panel legends removed',individual_figures_unchanged=True,scene_pixels_unchanged=True,scene_physical_sizes_unchanged=True,panels=panels,text_outside_page=0,text_bbox_overlaps=overlaps,geometry_camera_data_and_line_styles_unchanged=True,new_solver_runs=0,new_3D_renders=0)
    (OUT/('PREVIEW_LAYOUT_QA.json' if args.preview else 'LAYOUT_QA.json')).write_text(json.dumps(audit,indent=2))
    if not args.preview:
        (OUT/'README.txt').write_text('三图合并版v2：只保留居中一份共享图例。总标题28pt、视角标题24pt、坐标数字18pt、坐标轴标题20pt、图例20pt。PNG为600dpi，PDF/SVG文字可编辑。直接复用原PDF内无文字的三维场景，像素内容与原图一致，物理显示尺寸不变；只调整合并排版及文字。原三张单图和旧合并图均保留。\n')
        files=[p for p in OUT.iterdir() if p.suffix in ['.png','.pdf','.svg']];(OUT/'DELIVERY_SHA256.json').write_text(json.dumps({p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in files},indent=2))
    print(json.dumps(dict(status='PASS',fonts=SIZES,legends=1,overlaps=overlaps,output=str(OUT))),flush=True)

if __name__=='__main__':main()
