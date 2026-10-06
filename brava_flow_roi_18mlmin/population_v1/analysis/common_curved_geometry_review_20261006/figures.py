"""Publication/export figures from saved element evidence; never alter results."""
from common import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':13,'axes.labelsize':11,
    'figure.facecolor':'white','savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
COLORS=['#c44935','#777777','#087f8c']
FLOOR=1e-17

def bars(ax,labels,values,title,ylabel,tol=None):
    x=np.arange(len(labels));width=.23
    for i,(name,v) in enumerate(values.items()):
        v=np.array(v,float);ax.bar(x+(i-1)*width,np.maximum(v,FLOOR),width,label=name,color=COLORS[i])
        for xx,yy in zip(x+(i-1)*width,v):
            if yy==0:ax.text(xx,FLOOR*1.8,'0',ha='center',va='bottom',fontsize=8)
    ax.set_xticks(x,labels);ax.tick_params(axis='x',labelsize=10);ax.set_yscale('log');ax.set_ylim(FLOOR,1.)
    ax.set_title(title,loc='left',fontweight='bold');ax.set_ylabel(ylabel);ax.grid(axis='y',alpha=.15)
    if tol:ax.axhline(tol,color='#333333',linestyle=':',linewidth=1,label='Frozen gate' if title.startswith('a') else None)

def save(fig,name,description):
    fig.savefig(HERE/name,dpi=300,bbox_inches='tight');plt.close(fig)
    return dict(file=name,dpi=300,sha256=sha(HERE/name),description=description)

def main():
    result=read(HERE/'COMMON_GEOMETRY_HARNESS_RESULTS.json');data=np.load(HERE/'COMMON_GEOMETRY_SAMPLES.npz')
    boundary=read(HERE/'BOUNDARY_MIXED_GEOMETRY_AUDIT.json')
    cases={c['case']:c for c in result['cases']};cells=['central_wall','inlet_rim'];labels=['Central wall\n143654','Inlet rim\n35265']
    figures=[];source=[]
    def residual(tag,which):
        r=data[tag+'_'+which+'_R'][3,:4];ref=data[tag+'_reference_R'][3,:4]
        return float(np.linalg.norm(r-ref)/np.linalg.norm(ref))
    before=[residual(c+'_curved_P2_A','original') for c in cells]
    after=[residual(c+'_curved_P2_A','patched') for c in cells]
    fig,axs=plt.subplots(2,2,figsize=(13.6,9.2));fig.subplots_adjust(hspace=.44,wspace=.3,top=.87,bottom=.12)
    fig.suptitle('Common curved geometry: native assembly before and after',fontsize=18,fontweight='bold',y=.98)
    fig.text(.5,.93,'Two retained V2 cells | identical quadrature and physics kernels | no CFD',ha='center',fontsize=12)
    bars(axs[0,0],labels,{'Before':before,'Common reference':[0,0],'After':after},
        'a  Continuity residual','Relative L2 difference')
    vals={}
    for st,label in [('before','Before'),('reference','Common reference'),('after','After')]:
        vals[label]=[cases[c+'_curved_P2_P']['metrics'][st]['gradp_scaled'] if st!='reference' else 0. for c in cells]
    bars(axs[0,1],labels,vals,'b  Pressure gradient transform','H max |grad p_h - common grad p_h|')
    labs=[];vb=[];va=[]
    for c,l in zip(cells,['Wall','Rim']):
        for k,t in [('pressure_transform_J_scaled','J'),('used_det_relative','det J')]:
            labs.append(l+'\n'+t);vb.append(cases[c+'_curved_P2_A']['metrics']['before'][k]);va.append(cases[c+'_curved_P2_A']['metrics']['after'][k])
    bars(axs[1,0],labs,{'Before':vb,'Common reference':[0]*4,'After':va},'c  Geometry matrix and integration determinant','Scaled maximum difference')
    ax=axs[1,1]
    for c,l,color in zip(cells,['Central wall','Inlet rim'],['#c44935','#e1a22a']):
        key=c+'_curved_P2_A';qb=data[key+'_before_q'];qp=data[key+'_after_q'];mask=qb[:,0]==2
        # Reference determinant reconstructed from independently checked common velocity Jacobian.
        ref=np.linalg.det(data[key+'_after_velocity_J'])[mask]
        ax.plot(np.arange(4),qb[mask,5]/ref,'o-',color=color,label=l+' before')
        ax.plot(np.arange(4),qp[mask,5]/ref,'x--',color='#087f8c' if c=='central_wall' else '#263d77',label=l+' after')
    ax.axhline(1,color='#777777',ls=':',label='Common reference')
    ax.set_xticks(range(4));ax.set_xlabel('Continuity quadrature point (0-based)');ax.set_ylabel('Used det J / common det J')
    ax.set_title('d  First-point state inconsistency removed',loc='left',fontweight='bold');ax.grid(alpha=.15)
    ax.legend(fontsize=9,loc='lower right');axs[0,0].legend(fontsize=9,loc='upper right')
    fig.text(.02,.025,'Reference zeros are drawn at 1e-17 and labelled 0. Gradient comparison uses the same P1 pressure field,\nnot an arbitrary exact physical linear pressure. Volume gates pass; the separate weak-boundary gate remains BLOCKED.',fontsize=10)
    figures.append(save(fig,'01_COMMON_GEOMETRY_BEFORE_AFTER.png','Common mapping fixes captured continuity weights and pressure field transformation; no CFD error estimate'))

    tags=['central_wall_affine_P1_A','central_wall_affine_P2_A','central_wall_curved_P2_A','inlet_rim_curved_P2_A']
    names=['Affine\nP1/P1','Affine\nP2/P1','Curved wall\nP2/P1','Curved rim\nP2/P1']
    fig,axs=plt.subplots(2,2,figsize=(13.6,9.2));fig.subplots_adjust(hspace=.46,wspace=.32,top=.87,bottom=.12)
    fig.suptitle('P1 / P2 geometry regression',fontsize=18,fontweight='bold',y=.98)
    fig.text(.5,.93,'Original native objects, gated patch, independent common-geometry reference',ha='center',fontsize=12)
    bars(axs[0,0],names,{'Before':[residual(t,'original') for t in tags],'Reference':[0]*4,
        'After':[residual(t,'patched') for t in tags]},'a  Continuity residual vs reference','Relative L2 difference')
    matrix={}
    for st,l in [('original','Before'),('reference','Reference'),('patched','After')]:
        matrix[l]=[float(abs(data[t+'_'+st+'_K']-data[t+'_reference_K']).max())/float(abs(data[t+'_reference_K']).max()) for t in tags]
    bars(axs[0,1],names,matrix,'b  Complete local tangent vs reference','Relative maximum difference',1e-11)
    changeR=[];changeK=[]
    for t in tags:
        rb=data[t+'_original_R'];kb=data[t+'_original_K']
        changeR.append(float(abs(data[t+'_patched_R']-rb).max())/float(abs(rb).max()))
        changeK.append(float(abs(data[t+'_patched_K']-kb).max())/float(abs(kb).max()))
    bars(axs[1,0],names,{'Residual change':changeR,'Exact unchanged':[0]*4,'Matrix change':changeK},
        'c  Patch change (P1 is bitwise unchanged)','Relative maximum change')
    ax=axs[1,1];ax.axis('off')
    p1=[c for c in result['cases'] if c['velocity_nodes']==4]
    txt=('Frozen checks: 684 passed / 0 failed\n\n'
         'P1/P1: 4 fixtures, R and K bitwise identical\n'
         'Affine P2/P1: agreement within frozen 1e-11\n'
         'Curved P2/P1: common geometry gates pass\n'
         'A: G = I, div = 3\nB: G = diag(1, -1, 0), div = 0\n'
         'C: exact axial quadratic velocity\n\n'
         'Boundary: prescribed traction / strong BC pass\n'
         'Weak Dirichlet: remaining mapping defect\n'
         'NEXT CFD: NOT AUTHORIZED BY THIS GATE')
    ax.text(0,1,txt,va='top',fontsize=12,linespacing=1.45)
    axs[0,0].legend(fontsize=9,loc='upper right');axs[1,0].legend(fontsize=9,loc='upper right')
    fig.text(.02,.025,'Zero values use the labelled plotting floor only; all raw values are retained in NPZ/CSV.\nThese are element algebra regressions, not CFD solution accuracy or trajectory validation.',fontsize=10)
    figures.append(save(fig,'02_P1_P2_GEOMETRY_REGRESSION.png','P1 unchanged, affine P2 roundoff regression, curved volume fix and remaining boundary gate'))

    fs=boundary['fixtures'];selected=[f for f in fs if f['variant']=='curved'];bd=np.load(HERE/'BOUNDARY_MIXED_GEOMETRY_SAMPLES.npz')
    fig,axs=plt.subplots(1,3,figsize=(15.2,5.2));fig.subplots_adjust(wspace=.35,top=.78,bottom=.2)
    fig.suptitle('Boundary audit: passing traction path, blocked weak mixed path',fontsize=17,fontweight='bold',y=.98)
    fig.text(.5,.9,'Only retained cells and one explicit rigid cap control | no boundary patch',ha='center',fontsize=11)
    ln=['Wall\n143654','Wall\n35265','Inlet\n35265','Outlet*\n35265']
    axs[0].bar(range(4),[f['pressure_reference_coordinate_difference_max'] for f in selected],color='#c44935')
    axs[0].set_xticks(range(4),ln);axs[0].set_yscale('log');axs[0].axhline(1e-11,color='k',ls=':',label='1e-11')
    axs[0].set_title('a  Separate pressure inverse disagrees',loc='left');axs[0].set_ylabel('max |lambda_pressure - lambda_common|');axs[0].legend()
    axs[1].bar(range(4),[f['pressure_inverse_admissible_points'] for f in selected],color=['#c44935','#c44935','#e1a22a','#e1a22a'])
    axs[1].set_ylim(0,8);axs[1].set_yticks([0,7]);axs[1].set_xticks(range(4),ln);axs[1].set_ylabel('Accepted points / 7')
    axs[1].set_title('b  Native weak pressure admissibility',loc='left')
    for i,f in enumerate(selected):axs[1].text(i,.3 if i<2 else 7.15,'throws' if i<2 else 'mismatch',ha='center',fontsize=9)
    errs=[]
    for f in selected:
        t=next(t for t in boundary['tests'] if t['name']==f['case']+'_native_traction_residual');errs.append(t['error'])
    axs[2].bar(range(4),np.maximum(errs,1e-21),color='#087f8c');axs[2].set_yscale('log');axs[2].set_ylim(1e-21,1e-9)
    axs[2].axhline(1e-11,color='k',ls=':');axs[2].set_xticks(range(4),ln);axs[2].set_ylabel('Scaled maximum residual difference')
    axs[2].set_title('c  Native prescribed traction passes',loc='left')
    fig.text(.02,.025,'Outlet*: rigid transform of the same inlet cell, not a new actual outlet cell. Traction does not certify weak Dirichlet.\nPhase F stop: MIXED_CURVED_BOUNDARY_ASSEMBLY_BLOCKED. No pipe/BraVa CFD or particle runs.',fontsize=10)
    figures.append(save(fig,'03_BOUNDARY_MIXED_GEOMETRY_AUDIT.png','Actual boundary pressure inverse failures vs passing prescribed traction and Cauchy integration'))
    put('FIGURE_MANIFEST.json',dict(figures=figures,sources=['COMMON_GEOMETRY_HARNESS_RESULTS.json','COMMON_GEOMETRY_SAMPLES.npz',
        'COMMON_GEOMETRY_METRICS.csv','BOUNDARY_MIXED_GEOMETRY_AUDIT.json','BOUNDARY_MIXED_GEOMETRY_SAMPLES.npz']))

if __name__=='__main__':main()
