#!/usr/bin/env python3
"""Replot Particle-3 Figure 06 from saved data; preserve the historical figure."""
from pathlib import Path
import hashlib,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse,PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.text import Text
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/'particle_3d/reports/particle3'
SOURCE=REPORT/'data/06_deformation_examples.json'
ORIGINAL=REPORT/'figures/06_capillary_deformation_surrogate.png'
STEM='06_capillary_deformation_surrogate_academic_en'
BLUE='#0072B2';ORANGE='#D55E00';WALL='#64717C';INK='#202832'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def capsule_path(radius,length):
    upper=MplPath.arc(0,180);lower=MplPath.arc(180,360)
    u=upper.vertices*radius+[0,length/2]
    l=lower.vertices*radius-[0,length/2]
    codes=lower.codes.copy();codes[0]=MplPath.LINETO
    return MplPath(np.vstack([u,l,u[0]]),np.r_[upper.codes,codes,MplPath.CLOSEPOLY])


def main():
    before={str(p):digest(p) for p in [SOURCE,ORIGINAL]}
    source_receipt=json.loads((REPORT/'data/06_figure_sources.json').read_text())
    assert before[str(SOURCE)]==source_receipt['source_sha256'][SOURCE.name]
    rows=json.loads(SOURCE.read_text());rid=rows[0]['rbc_id']
    same=[r for r in rows if r['rbc_id']==rid]
    statuses=['FREE_OBLATE','CAPILLARY_DEFORMED','DEFORMATION_SURROGATE_INFEASIBLE']
    selected=[next(r for r in same if r['status']==s) for s in statuses]
    assert all(r['V_fL']==selected[0]['V_fL'] and r['area_budget_m2']==selected[0]['area_budget_m2'] for r in selected)
    deformed=selected[1];r=deformed['R_cap_m'];length=deformed['L_cap_m']
    volume=np.pi*r*r*length+4*np.pi*r**3/3
    area=2*np.pi*r*length+4*np.pi*r*r
    assert np.isclose(volume,deformed['volume_m3'],rtol=1e-12,atol=0)
    assert np.isclose(area/deformed['area_budget_m2'],deformed['area_ratio'],rtol=1e-12,atol=0)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'mathtext.fontset':'dejavusans',
        'axes.labelsize':9.5,'axes.linewidth':.65,'axes.edgecolor':INK,'text.color':INK,
        'xtick.labelsize':8.5,'ytick.labelsize':8.5,'xtick.color':INK,'ytick.color':INK,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
        'svg.fonttype':'none','savefig.facecolor':'white'})
    fig=plt.figure(figsize=(7.6,6.4),facecolor='white')
    axes=fig.subplots(1,3,sharex=True,sharey=True)
    fig.subplots_adjust(left=.085,right=.99,bottom=.322,top=.796,wspace=.12)
    fig.text(.5,.975,'Geometric feasibility of an RBC deformation surrogate',
        ha='center',va='top',fontsize=12.5,weight='bold')
    fig.text(.5,.928,
        rf"Same reference cell: $V_0 = {selected[0]['V_fL']:.3f}$ fL; "
        rf"$A_{{\mathrm{{max}}}} = {selected[0]['area_budget_m2']*1e12:.3f}$ µm²",
        ha='center',va='center',fontsize=9.2)
    handles=[Line2D([],[],color=BLUE,lw=1.3,ls=(0,(4,2)),label='Original oblate geometry'),
        Line2D([],[],color=ORANGE,lw=1.5,label='Selected capsule'),
        Line2D([],[],color=WALL,lw=1.5,label='Tube wall')]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,.887),ncol=3,
        frameon=False,fontsize=8.7,handlelength=2.4,columnspacing=2.)
    titles=['Free oblate','Feasible capsule','Surrogate infeasible']
    for i,(ax,record,title) in enumerate(zip(axes,selected,titles)):
        halfwidth=record['tube_apothem_m']*1e6
        ax.axvspan(-halfwidth,halfwidth,color='#F3F5F7',lw=0,zorder=0)
        for sign in [-1,1]:ax.axvline(sign*halfwidth,color=WALL,lw=1.35,zorder=2)
        if 'R_cap_m' in record:
            path=capsule_path(record['R_cap_m']*1e6,record['L_cap_m']*1e6)
            ax.add_patch(PathPatch(path,facecolor=ORANGE,alpha=.10,edgecolor='none',zorder=3))
            ax.add_patch(PathPatch(path,facecolor='none',edgecolor=ORANGE,lw=1.5,zorder=4))
        ax.add_patch(Ellipse((0,0),2*record['a_m']*1e6,2*record['c_m']*1e6,
            facecolor='none',edgecolor=BLUE,ls=(0,(4,2)),lw=1.25,zorder=5))
        ax.set(xlim=(-5.2,5.2),ylim=(-7.5,7.5),xticks=[-4,-2,0,2,4],yticks=[-6,-3,0,3,6])
        ax.set_aspect('equal',adjustable='box');ax.tick_params(direction='out',length=3,width=.65,pad=3)
        ax.text(.5,1.110,f'({chr(97+i)})  {title}',transform=ax.transAxes,
            ha='center',va='center',fontsize=9.6,weight='bold')
        ax.text(.5,1.037,rf'Wall half-width $h = {halfwidth:.3f}$ µm',transform=ax.transAxes,
            ha='center',va='center',fontsize=8.6,color='#4E5964')
    axes[0].set_ylabel('Axial position, z (µm)',labelpad=6)
    fig.text(.535,.267,'Transverse position, x (µm)',ha='center',fontsize=9.5)
    tabax=fig.add_axes([.085,.064,.905,.164]);tabax.set_axis_off()
    body=[['Wall clearance (µm)',f"{selected[0]['gap_m']*1e6:.3f}",f"{deformed['gap_m']*1e6:.3f}",'—'],
        ['Capsule radius, R (µm)','—',f"{deformed['R_cap_m']*1e6:.3f}",'—'],
        ['Cylinder length, L (µm)','—',f"{deformed['L_cap_m']*1e6:.3f}",'—'],
        ['Capsule area / area budget','—',f"{deformed['area_ratio']:.3f}",'—']]
    table=tabax.table(cellText=body,colLabels=['Geometric quantity','(a) Free','(b) Feasible','(c) Infeasible'],
        cellLoc='center',colLoc='center',colWidths=[.40,.20,.20,.20],bbox=[0,0,1,1])
    table.auto_set_font_size(False);table.set_fontsize(8.3)
    for (row,col),cell in table.get_celld().items():
        cell.visible_edges='';cell.set_facecolor('white');cell.PAD=.025
        if col==0:cell.get_text().set_ha('left')
        if row==0:
            cell.visible_edges='TB';cell.set_edgecolor('#939DA5');cell.set_linewidth(.65)
            cell.get_text().set_fontweight('bold')
        if row==len(body):cell.visible_edges='B';cell.set_edgecolor('#C2C8CD');cell.set_linewidth(.5)
    fig.text(.5,.023,'Geometric surrogate only; infeasibility does not imply physiological occlusion.',
        ha='center',va='center',fontsize=8.1,color='#4E5964')
    fig.canvas.draw();renderer=fig.canvas.get_renderer();extent=fig.bbox
    labels=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    assert not any(re.search(r'[\u3400-\u9fff]',t.get_text()) for t in labels)
    overflow=[]
    for t in labels:
        b=t.get_window_extent(renderer)
        if b.x0<extent.x0-1 or b.y0<extent.y0-1 or b.x1>extent.x1+1 or b.y1>extent.y1+1:overflow.append(t.get_text())
    assert not overflow,overflow
    scales=[]
    for ax in axes:
        p=ax.transData.transform([[0,0],[1,0],[0,1]])
        scales.append([float(np.linalg.norm(p[1]-p[0])),float(np.linalg.norm(p[2]-p[0]))])
    assert np.ptp(np.array(scales))<1e-9
    outputs=[]
    for fmt in ['png','pdf','svg']:
        p=REPORT/'figures'/f'{STEM}.{fmt}'
        kwargs={'dpi':600} if fmt=='png' else {}
        fig.savefig(p,**kwargs);outputs.append(dict(path=str(p),sha256=digest(p)))
    plt.close(fig)
    with Image.open(outputs[0]['path']) as im:dimensions=list(im.size);im.verify()
    assert dimensions==[4560,3840]
    assert all(digest(p)==value for p,value in before.items())
    caption=f'''Geometric feasibility of a volume-constrained red blood cell (RBC) deformation surrogate.
The same original oblate cell is used in all panels (cell ID {rid}; volume {selected[0]['V_fL']:.3f} fL;
surface-area budget {selected[0]['area_budget_m2']*1e12:.3f} µm²). (a) The original cell is geometrically admissible
in the wider tube. (b) In the narrower tube, a selected capsule preserves cell volume and satisfies
the surface-area inequality A <= A_max. L denotes the cylindrical segment length, excluding the
two hemispherical end caps. (c) No capsule satisfies the surrogate constraints in the narrowest tube.
Blue dashed lines denote the original oblate geometry; orange outlines denote the selected capsule.
Grey lines show walls at +/-h, using the polygonal tube apothem as the schematic wall half-width;
the nominal circumradius is not substituted for h. All panels have identical axis limits and equal
geometric scaling. The dash in the table means that the selected-shape quantity is not applicable
or no admissible shape was selected. The model imposes a surface-area budget, not surface-area
equality, and does not solve membrane mechanics. Infeasibility is a geometric surrogate outcome,
not evidence of physiological occlusion. These artificial tubes are validation geometries, not
a physiological capillary-size distribution.
'''
    (REPORT/f'{STEM}_caption.txt').write_text(caption)
    receipt=dict(all_pass=True,source_data=str(SOURCE),source_sha256=before[str(SOURCE)],
        original_figure=str(ORIGINAL),original_figure_sha256=before[str(ORIGINAL)],
        historical_source_and_figure_unchanged=True,selected_rbc_id=rid,
        selection='SAME_RECORDS_AS_ORIGINAL_STAGE6_FIRST_MATCH_PER_STATUS',selected_records=selected,
        changes=['English-only labels','Shared equal geometric scale','Shared legend and reference parameters',
            'Aligned numerical table outside the geometry','Tube apothem explicitly distinguished from nominal radius',
            '600 dpi PNG and editable vector PDF/SVG'],
        image_dimensions=dimensions,png_dpi=600,figure_size_inches=[7.6,6.4],
        pixel_scale_at_layout_dpi=scales,no_text_outside_canvas=True,english_labels_checked=True,
        capsule_volume_relative_error=abs(volume-deformed['volume_m3'])/deformed['volume_m3'],
        capsule_area_ratio=area/deformed['area_budget_m2'],outputs=outputs,
        script_sha256=digest(__file__),manual_visual_review='PENDING_USER_REVIEW')
    (REPORT/'data'/f'{STEM}_validation.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(all_pass=True,dimensions=dimensions,outputs=outputs),indent=2))


if __name__=='__main__':main()
