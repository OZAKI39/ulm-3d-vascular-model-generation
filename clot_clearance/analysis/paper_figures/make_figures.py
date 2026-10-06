"""Render archived scientific evidence. Does not import or call any solver."""
import argparse
import csv
import hashlib
import json
import platform
import re
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch
from matplotlib.text import Text

import common_style as style
from common_style import plt, COLORS as C, MARKERS as M, MM, panels, legend
import data_loading as d
from figure_specs import SPECS, write_plan, V, R, B, A, C5, C7, ONSET, S5, S7, BS

HERE=d.HERE
OUT=HERE/'output'
FIG=OUT/'figures'
MAN=OUT/'manifests'
NARRATIVE=[]

def curve(ax,x,y,c,*,label=None,open_marker=False,**kwargs):
    ax.plot(x,y,color=C[c],marker=M[c],markerfacecolor='white' if open_marker else C[c],
            label=label or c,**kwargs)

def note(ax,text,y=.97,**kwargs):
    ax.text(.03,y,text,transform=ax.transAxes,va='top',fontsize=7,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=2),**kwargs)

def title(fig,text):
    fig.suptitle(text,fontsize=10)

def record_table(name,rows):
    d.table(name+'.csv',rows)

def overview():
    ons=d.onset(); br=d.branches(); b=d.boundary(); status=d.read_json(V+B+'SUMMARY.json')
    f,ax=plt.subplots(figsize=(180*MM,121*MM));ax.set_axis_off();ax.set(xlim=(0,1),ylim=(0,1))
    title(f,'From bond activation to interpretable clot mechanics')
    positive={c:[r for r in br if r['configuration']==c and r['positive']] for c in ('C5','C7')}
    rows=[
        ('(a) Observed activation',C['C5'],
         f"C5: {ons[0]['strain']*100:.4f}% strain; {ons[0]['active_count']} / {ons[0]['total_bonds']:,} bonds\n"
         f"C7: {ons[1]['strain']*100:.4f}% strain; {ons[1]['active_count']} / {ons[1]['total_bonds']:,} bonds",
         'Sampling changes the apparent fraction.\nMesh objectivity is not established.'),
        ('(b) Numerical evidence',C['C7'],
         f"C5: {len(positive['C5'])} audited stable positive states\n"
         f"C7: {len(positive['C7'])} stationary positive states; stability unassessed",
         f"Half-reaction / continuum: {b['C5']['top_required']['ratio_to_continuum']:.3f} (C5),\n"
         f"{b['C7']['top_required']['ratio_to_continuum']:.3f} (C7); loading issue unresolved."),
        ('(c) Clot-model implications',C['L050'],
         'Validate finite-domain loading and force support.\nDefine a physically interpretable failure observable.',
         'Then test clot-specific calibration and coupling.\nCurrent coupon evidence is not a rupture prediction.')]
    for i,(head,col,left,right) in enumerate(rows):
        y=.69-i*.30
        ax.add_patch(FancyBboxPatch((.02,y),.96,.235,boxstyle='round,pad=0.006,rounding_size=.008',
                                  facecolor='#f7f8fa',edgecolor='#c7cbd1',linewidth=.7))
        ax.text(.04,y+.19,head,color=col,fontsize=10,weight='bold',va='center')
        ax.text(.04,y+.118,left,fontsize=8,va='top',linespacing=1.65)
        ax.text(.555,y+.118,right,fontsize=8,va='top',linespacing=1.65)
        ax.plot([.53,.53],[y+.035,y+.135],color='#d0d4da',lw=.8)
        if i<2:ax.annotate('',xy=(.5,y-.053),xytext=(.5,y-.008),arrowprops=dict(arrowstyle='->',color='#888888'))
    ax.text(.02,.005,'Research argument only; arrows do not establish causality or a biological mechanism.',fontsize=7,color=C['gray'])
    f.subplots_adjust(left=.015,right=.985,bottom=.035,top=.91)
    assert status['status']['ROOT_CAUSE_PRIMARY']=='INCONCLUSIVE'
    return f

def activation():
    rows=d.onset();f,ax=panels(2,2,height=120)
    title(f,'Activation onset and bond sampling')
    for i,r in enumerate(rows):
        c=r['configuration'];v=100*r['strain']
        ax[0].errorbar(i,v,yerr=[[100*(r['strain']-r['bracket_lower'])],[100*(r['bracket_upper']-r['strain'])]],
                      fmt=M[c],color=C[c],capsize=5,ms=5)
        ax[0].annotate(f'{v:.4f}%',(i,v),xytext=(0,9),textcoords='offset points',ha='center',fontsize=8)
        v=100*r['active_fraction'];ax[1].bar(i,v,color=C[c],width=.5)
        ax[1].text(i,v+.002,f"{r['active_count']} / {r['total_bonds']:,}\n{v:.5f}%",ha='center',va='bottom',fontsize=7)
        for j,(key,lab) in enumerate([('h_m','Spacing'),('median_depth_m','Activation depth')]):
            ax[2].bar(i+(j-.5)*.28,r[key]*1e6,width=.26,color=C[c],alpha=1 if j==0 else .4,
                      hatch=None if j==0 else '///',label=lab if i==0 else None)
        ax[3].bar(i,r['depth_over_delta'],width=.5,color=C[c])
        ax[3].text(i,r['depth_over_delta']+.007,f"{r['depth_over_delta']:.4f}\n$d/h={r['depth_over_h']:.1f}$",ha='center',fontsize=7)
    for a in ax:a.set_xticks([0,1],['C5','C7']);a.set_xlim(-.55,1.55)
    ax[0].set(ylabel='Diagnostic onset strain (%)',ylim=(.28,.49),title='Onset strain')
    note(ax[0],'Intervals: numerical brackets',y=.14)
    ax[1].set(ylabel='First-set bond fraction (%)',ylim=(0,.049),title='Same count, different denominator')
    ax[2].set(ylabel='Length ($\\mu$m)',ylim=(0,85),title='Sampling and activation location');legend(ax[2])
    ax[3].set(ylabel=r'Median clamp distance $d/\delta$',ylim=(0,.165),title='Activation remains near the clamps')
    return f

def precursor():
    b=d.boundary();comp=d.clamp();config=d.read_json(V+A+'CONFIG.json')
    eps=config['phase_A_strains'][-1];assert eps==.003
    f,ax=panels(2,2,height=127);title(f,f'Intact equilibria at engineering strain {eps:g}')
    gaprows=[];driverrows=[];control=[]
    for c in ('C5','C7'):
        rows=d.gaps(c);z=(d.numeric(rows,'z_left_m')+d.numeric(rows,'z_right_m'))/(2*b[c]['L_m'])
        gap=d.numeric(rows,'displacement_gap_m')*1e6;strain=d.numeric(rows,'engineering_strain')*100
        assert np.allclose(d.numeric(rows,'displacement_gap_m')/d.numeric(rows,'reference_gap_m'),strain/100)
        curve(ax[0],z,gap,c);curve(ax[1],z,strain,c)
        for r,zz,gg,ss in zip(rows,z,gap,strain):
            gaprows.append(dict(configuration=c,z_mid_over_L=zz,displacement_increment_um=gg,local_engineering_strain_percent=ss,**r))
        bonds=d.read_csv(R+A+f'tables/{c}_FULL_CLAMP_BASELINE_0.0030_AXIAL_BONDS.csv')
        xx=d.numeric(bonds,'d_over_delta');yy=d.numeric(bonds,'A_psi')
        ax[3].scatter(xx,yy,s=12,color=C[c],marker=M[c],facecolors='none',alpha=.55,label=c,linewidths=.6)
        for r in bonds:driverrows.append(dict(configuration=c,**r))
        x=next(x for x in comp if x['level']==int(c[1:]));assert np.isclose(yy.max(),x['baseline_max_Apsi'])
        control.append(dict(configuration=c,baseline_max_Apsi=x['baseline_max_Apsi'],axial_only_max_Apsi=x['control_max_Apsi']))
    for j,(key,lab,hatch) in enumerate([('baseline_max_Apsi','Full clamp',None),('axial_only_max_Apsi','Axial only','///')]):
        v=[r[key] for r in control];bars=ax[2].bar(np.arange(2)+(j-.5)*.30,v,.28,color=[C['C5'],C['C7']],alpha=1 if j==0 else .45,hatch=hatch,label=lab)
        ax[2].bar_label(bars,fmt='%.2f',padding=3,fontsize=7)
    ax[0].set(xlabel=r'Centerline midpoint $z/L$',ylabel=r'Displacement increment ($\mu$m)',title='Centerline displacement increments');legend(ax[0],loc='lower center')
    ax[1].axhline(eps*100,color=C['gray'],ls='--',label='Affine reference');legend(ax[1])
    ax[1].set(xlabel=r'Centerline midpoint $z/L$',ylabel='Local engineering strain (%)',title='Local strain increases as spacing falls')
    ax[2].set(xticks=[0,1],xticklabels=['C5','C7'],ylim=(0,14.8),ylabel=r'Maximum $A_\psi=\Psi^+_{eq}/\Psi^+_{aff}$',title='Lateral release does not reduce the peak');legend(ax[2])
    ax[3].set(xlabel=r'Bond distance to clamp $d/\delta$',ylabel=r'$A_\psi=\Psi^+_{eq}/\Psi^+_{aff}$',title='Saved axial-bond amplification');legend(ax[3])
    record_table('fig03_centerline',gaprows);record_table('fig03_axial_bonds',driverrows);record_table('fig03_loading_control',control)
    return f

def layer_forces():
    b=d.boundary();f,ax=panels(1,2,height=90);title(f,'Affine forces and manufactured axial loading')
    layers=[];cum=[]
    for c in ('C5','C7'):
        rows=d.layer(c);x=d.numeric(rows,'z_m')/b[c]['L_m'];y=d.numeric(rows,'layer_fz_over_R')
        assert np.allclose(y,d.numeric(rows,'sum_fz_N')/b[c]['continuum']['reaction_N'])
        curve(ax[0],x,y,c)
        layers.extend([dict(configuration=c,z_over_L=xx,**r) for xx,r in zip(x,rows)])
    for c in ('C5','C7','L075'):
        rows=d.cumulative(c);curve(ax[1],d.numeric(rows,'depth_over_delta'),d.numeric(rows,'cumulative_over_R'),c)
        cum.extend([dict(configuration=c,**r) for r in rows])
    ax[0].axhline(0,color=C['gray'],lw=.65)
    ax[0].set(xlabel=r'Layer position $z/L$',ylabel=r'Internal layer force $\sum_i f_{iz}/R_{cont}$',title='Signed layer forces');legend(ax[0])
    ax[1].axhline(1,color=C['gray'],ls='--',label='Continuum reaction')
    ax[1].set(xlabel=r'Depth from top surface $d/\delta$',ylabel=r'Cumulative required force / $R_{cont}$',title='Top-half cumulative external force',ylim=(.25,1.65));legend(ax[1],loc='lower right')
    note(ax[1],'Manufactured load: $f^{ext}=-f^{int}$\n$R_{cont}=P_{zz}L^2$',y=.97)
    record_table('fig04_layers',layers);record_table('fig04_cumulative',cum)
    return f

def components():
    b=d.boundary();norm=d.read_json(V+B+'NORMALIZATION_CHANNELS.json')['cases'];f,ax=panels(2,2,height=126)
    title(f,'First-interior center-particle force: exact signed decompositions')
    output=[]
    pairs=[('Force formula',['Direct','Z correction']),('Neighbor location',['Loaded surface','Other neighbors']),
           ('Endpoint normalization',['Loaded surface','Other endpoints']),('Normalization channels',['Lattice-weight\ncomponent','Normalization\nexcess'])]
    for i,c in enumerate(('C5','C7')):
        r=next(r for r in b[c]['first_layer'] if r['side']=='bottom')
        n=next(x for x in norm if x['case']==c);n=next(x for x in n['rows'] if x['ID']==r['particle_ID'])
        vals=[[r['direct_fraction'],r['Z_fraction']],[r['boundary_neighbor_fraction'],r['interior_neighbor_fraction']],
              [r['endpoint_signed_fractions']['loaded_surface_endpoint'],r['endpoint_signed_fractions']['nonloaded_endpoint']],
              [n['uniform_wbar_component_signed_fraction'],n['normalization_excess_component_signed_fraction']]]
        for a,(head,labs),v in zip(ax,pairs,vals):
            assert abs(sum(v)-1)<1e-12
            bars=a.bar(np.arange(2)+(i-.5)*.32,v,.30,color=C[c],label=c)
            a.bar_label(bars,labels=[f'{x:+.3f}' for x in v],padding=3,fontsize=7)
            a.set(xticks=[0,1],xticklabels=labs,ylabel=r'Component force / net $f_z$',title=head,ylim=(-.39,1.48))
            output.extend([dict(configuration=c,particle_ID=r['particle_ID'],decomposition=head,component=lab.replace('\n',' '),signed_fraction=vj) for lab,vj in zip(labs,v)])
    for a in ax:
        a.axhline(0,color=C['gray'],lw=.65);a.axhline(1,color=C['gray'],ls=':',lw=.7)
    legend(ax[0],loc='upper left');record_table('fig05_signed_components',output)
    return f

def support():
    b=d.boundary();cfg=d.read_json(V+B+'CONFIG.json');f,ax=panels(2,2,height=129)
    title(f,'Geometric family support and affine force residuals')
    rowsout=[]
    for a,c in zip(ax,('C5','L050','L075')):
        rows=d.particles(c);x=d.numeric(rows,'d_over_delta');y=d.numeric(rows,'force_norm_over_R');full=np.array([r['full_support']=='True' for r in rows])
        assert full.sum()==b[c]['full_support_count'];assert np.all(y>0)
        a.scatter(x[~full],y[~full],s=13,facecolors='none',edgecolors=C['gray'],alpha=.55,linewidths=.6,label='Truncated family')
        if full.any():a.scatter(x[full],y[full],s=21,color=C[c],marker='D',label='Full family',zorder=3)
        a.axhline(cfg['bulk_noise_normalized_tolerance'],color=C['C7'],ls='--',lw=.8,label='Core noise criterion')
        a.axvline(1,color=C['gray'],ls=':',lw=.8)
        a.set(xlabel=r'Nearest-surface distance $d_{min}/\delta$',ylabel=r'Nodal force $|f_i|/R_{cont}$',yscale='log',ylim=(1e-18,.3),xlim=(-.04,1.57),title=f"{c}: $L={b[c]['L_m']*1e3:g}$ mm")
        a.set_yticks([1e-16,1e-12,1e-8,1e-4]);legend(a,loc='lower left',fontsize=6.5)
        rowsout.extend([dict(configuration=c,**r) for r in rows])
    a=ax[3];a.axis('off')
    a.set_title('Available core and declared gates',pad=9)
    a.text(.02,.94,'Case        1-hop family     2-hop     3-hop',fontsize=8,transform=a.transAxes,va='top',family='DejaVu Sans Mono')
    for i,c in enumerate(('C5','C7','L050','L075')):
        q=b[c];a.text(.02,.82-i*.115,f"{c:<5} {q['full_support_count']:>12} {q['two_hop_closed_count']:>9} {q['three_hop_closed_count']:>9}",fontsize=8,transform=a.transAxes,family='DejaVu Sans Mono')
    q=b['L075'];a.text(.02,.32,
        f"L075 full-family core: {q['bulk_patch_test']}\n"
        f"max |f| / R = {q['bulk']['max_over_R']:.2e}  (limit {cfg['bulk_noise_normalized_tolerance']:.0e})\n"
        f"core / collar = {q['bulk_collar_max_ratio']:.4f}  (limit {cfg['bulk_collar_ratio_max']:.0e})\n"
        'Geometric completeness does not close\nthe complete force-dependency stencil.',
        transform=a.transAxes,fontsize=8,va='top',linespacing=1.6)
    record_table('fig06_particles',rowsout)
    record_table('fig06_support_summary',[dict(configuration=c,particles=r['particles'],full_family=r['full_support_count'],two_hop=r['two_hop_closed_count'],three_hop=r['three_hop_closed_count'],core_max_over_R=r['bulk']['max_over_R'],core_collar_ratio=r['bulk_collar_max_ratio'],noise_criterion=cfg['bulk_noise_normalized_tolerance'],ratio_criterion=cfg['bulk_collar_ratio_max'],patch_status=r['bulk_patch_test']) for c,r in b.items()])
    return f

def work_energy():
    b=d.boundary();cfg=d.read_json(V+B+'CONFIG.json');work=d.read_json(V+B+'VIRTUAL_WORK_TRACTION.json')
    f,ax=panels(1,3,height=89);title(f,'Global energy, affine virtual work and half-reaction')
    rows=[]
    for i,(c,r) in enumerate(b.items()):
        ee=abs(r['global_affine_energy_ratio']-1);we=abs(r['affine_virtual_work_over_RL']-1);rr=r['top_required']['ratio_to_continuum']
        assert ee>0 and we>0
        for a,v in zip(ax,(ee,we,rr)):a.plot(i,v,marker=M[c],color=C[c],ms=5)
        ax[2].annotate(f'{rr:.3f}',(i,rr),xytext=(0,7),textcoords='offset points',ha='center',fontsize=7)
        ww=next(x for x in work if x['case']==c);assert np.isclose(rr,ww['half_force_ratio'])
        rows.append(dict(configuration=c,energy_ratio=r['global_affine_energy_ratio'],absolute_relative_energy_error=ee,affine_virtual_work_ratio=r['affine_virtual_work_over_RL'],absolute_relative_work_error=we,half_reaction_ratio=rr,force_centroid_m=ww['signed_force_centroid_z_m'],settling_depth_m=ww['final_signed_cumulative_settled_1pct_depth_m']))
    for a in ax:a.set_xticks(range(4),list(b),rotation=30);a.set_xlim(-.5,3.5)
    ax[0].set(yscale='log',ylim=(1e-15,1e-12),ylabel=r'$|U_{PD}/U_{cont}-1|$',title='Affine energy error')
    ax[1].set(yscale='log',ylim=(1e-15,1e-12),ylabel=r'$|\mathcal{W}_{aff}/(R_{cont} L)-1|$',title='Virtual-work error')
    tol=cfg['traction_relative_tolerance'];ax[2].axhspan(1-tol,1+tol,color='#dddddd',label=f'{tol*100:g}% audit band')
    ax[2].axhline(1,color=C['gray'],lw=.7,ls='--');ax[2].set(ylim=(.93,1.57),ylabel=r'$R_{half}/R_{cont}$',title='Required half-reaction');legend(ax[2],loc='upper left')
    record_table('fig07_energy_work_reaction',rows)
    return f

def positive_branch():
    onset={r['configuration']:r['strain'] for r in d.onset()};rows=d.branches();f,ax=panels(2,2,height=126)
    title(f,'Archived positive-phase branches')
    exported=[]
    for c in ('C5','C7'):
        rr=[r for r in rows if r['configuration']==c and r['positive']];x=d.numeric(rr,'strain');relative=(x-onset[c])*1e4
        for a,key in zip(ax,['maximum_s','weighted_mean_s','nominal_stress_Pa','active_count']):
            xx=x*100 if key=='nominal_stress_Pa' else relative
            curve(a,xx,d.numeric(rr,key),c,open_marker=c=='C7',label=c+(' (audited stable)' if c=='C5' else ' (stability unassessed)'))
        exported.extend([dict(strain_above_onset_times_1e4=float(dx),**r) for dx,r in zip(relative,rr)])
    for a in (ax[0],ax[1],ax[3]):a.set_xlabel(r'Strain above diagnostic onset $(\epsilon-\epsilon_{act})\times10^4$')
    ax[0].set(ylabel=r'Maximum phase $s$',title='Local phase amplitude');legend(ax[0],loc='upper left',fontsize=6.5)
    ax[1].set(ylabel=r'Quadrature-weighted mean phase $\bar{s}$',title='Domain-weighted phase');ax[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax[2].set(xlabel='Engineering strain (%)',ylabel='Nominal reaction stress (Pa)',title='Unequal load windows');legend(ax[2],loc='center left',fontsize=6.5)
    ax[3].set(ylabel='Positive-phase bond count',ylim=(0,3),yticks=[0,1,2,3],title='Support count remains localized')
    note(ax[3],'C5: filled markers; C7: open markers\nOverlapping count curves are both retained.',y=.3)
    record_table('fig08_positive_states',exported)
    return f

def stability():
    rows=d.branches();states=d.read_json(S5)['states'];tol=d.read_json(V+'pfpd_at1_petsc_vi_solver_v1/provenance/PREDECLARED_CONFIG.json');c7tol=d.read_json(V+C7+'provenance/PREDECLARED_CONFIG.json')
    cfg=d.read_json(V+C5+'provenance/PREDECLARED_CONFIG.json')
    assert c7tol['force_tolerance']==tol['force_tolerance'] and c7tol['KKT_tolerance']==tol['KKT_tolerance']
    f,ax=panels(2,2,height=129);title(f,'Stationarity and stability of the audited states')
    x=np.array([s['strain']*100 for s in states]);h=np.array([s['lambda_x_min'] for s in states]);sel=[s for s in states if s['lambda_S_min'] is not None]
    ax[0].plot(x,h,'o-',color=C['C5'],label=r'Mechanical $\lambda_{min}(H_{xx})$')
    ax[0].plot([s['strain']*100 for s in sel],[s['lambda_S_min'] for s in sel],'s-',color=C['L050'],label=r'Condensed $\lambda_{min}(S)$')
    thresholds=np.array([max(v['curvature_tol'] for v in s['eigenvalues']) for s in states])
    ax[0].plot(x,thresholds,'--',color=C['gray'],label='Archived curvature threshold')
    ax[0].set(yscale='log',ylim=(2e-8,.3),xlabel='Engineering strain (%)',ylabel='Dimensionless curvature',title='C5 physical-energy curvature');legend(ax[0],loc='center right',fontsize=6.5)
    note(ax[0],'Intact state: no movable phase block',y=.10)
    export=[]
    for c in ('C5','C7'):
        rr=[r for r in rows if r['configuration']==c];xx=d.numeric(rr,'index')+1
        for a,key,tkey in [(ax[1],'force_residual','force_tolerance'),(ax[2],'projected_phase_KKT','KKT_tolerance')]:
            y=d.numeric(rr,key)/tol[tkey];valid=y>0
            curve(a,xx[valid],y[valid],c,open_marker=c=='C7')
        export.extend([dict(force_over_tolerance=r['force_residual']/tol['force_tolerance'],KKT_over_tolerance=r['projected_phase_KKT']/tol['KKT_tolerance'],**r) for r in rr])
    for a in ax[1:3]:
        a.axhline(1,color=C['gray'],ls='--',label='Acceptance limit');a.set(yscale='log',xlabel='Saved-state index within each grid',ylim=(1e-11,4),xticks=[1,5,10,15,20]);legend(a,loc='lower right')
    ax[1].set(ylabel=r'$r_f/\tau_f$  ($\tau_f=10^{-9}$)',title='Force stationarity')
    ax[2].set(ylabel=r'$r_{KKT}/\tau_{KKT}$  ($\tau_{KKT}=10^{-8}$)',title='Projected phase stationarity')
    zero=sum(r['projected_phase_KKT']==0 for r in rows);note(ax[2],f'{zero} exact-zero value omitted on log axis',y=.93)
    npos=sum(s['metrics']['maximum_s']>cfg['s_tol'] for s in states)
    ax[3].axis('off');ax[3].set_title('What the evidence establishes',pad=10)
    ax[3].text(.03,.91,f"C5: {len(states)} / {len(states)} audited states stable\n"
               f"Includes {npos} positive-phase states.\nNo instability identified in this window.\n\n"
               f"C7: {sum(r['configuration']=='C7' for r in rows)} stationary states\nStability: NOT ASSESSED\nUnclassified does not mean unstable.\n\n"
               'No instability interval or rupture\ncriterion has been established.',transform=ax[3].transAxes,va='top',fontsize=8.5,linespacing=1.65)
    record_table('fig09_stationarity',export)
    record_table('fig09_C5_curvatures',[dict(label=s['label'],strain=s['strain'],lambda_Hxx_min=s['lambda_x_min'],lambda_S_min=s['lambda_S_min'],curvature_threshold=threshold,status=s['status'],movable_phase_count=s['m']) for s,threshold in zip(states,thresholds)])
    return f

def implications():
    ons=d.onset();br=d.branches();b=d.boundary();status=d.read_json(V+B+'SUMMARY.json');d.read_json(V+A+'SUMMARY.json')
    f,ax=plt.subplots(figsize=(180*MM,133*MM));ax.axis('off');title(f,'What current evidence implies for future clot simulation')
    stable=sum(r['configuration']=='C5' and r['positive'] for r in br)
    content=[
        [f"Two first-set bonds on each grid;\nonset shifts by {(ons[1]['strain']/ons[0]['strain']-1)*100:.1f}% (C7 vs C5).",'Bond counts and sampling\nare not a physical crack area.','Define and validate a failure\nobservable before rupture claims.'],
        [f"C5: {stable} stable positive states;\nC7 stability remains unassessed.",'Phase growth does not establish\ninstability or propagation.','Keep stationarity, stability and\nphysical failure tests separate.'],
        [f"Half-reaction ratios:\nC5 {b['C5']['top_required']['ratio_to_continuum']:.3f}; C7 {b['C7']['top_required']['ratio_to_continuum']:.3f}.",'Global energy agreement does\nnot validate physical traction.','Validate finite-domain loading\nand its work-conjugate forces.'],
        [f"L075: {b['L075']['full_support_count']} full-family particles;\n{b['L075']['three_hop_closed_count']} three-hop closed particles.",'The geometric core test fails;\nroot cause is inconclusive.','Use a complete force-support\noracle before assigning a cause.'],
        ['Current results are numerical\ncoupon and audit evidence.','No clot-specific calibration or\nbiological rupture is established.','Relate validated mechanics to\nclot data, then assess coupling.']]
    widths=[.34,.33,.33]
    tab=ax.table(cellText=content,colLabels=['Existing evidence','Interpretation limit','Development implication'],
                 cellLoc='left',colLoc='left',colWidths=widths,bbox=[.0,.08,1,.88])
    tab.auto_set_font_size(False);tab.set_fontsize(8)
    for (row,col),cell in tab.get_celld().items():
        cell.set_edgecolor('#d5d9dd');cell.set_linewidth(.65);cell.PAD=.065
        if row==0:cell.set_facecolor('#e6ecf2');cell.get_text().set_weight('bold')
        else:cell.set_facecolor('#fafbfc' if row%2 else 'white')
    ax.text(.005,.025,'Development implications, not completed experiments or patient-level predictions.',fontsize=7,color=C['gray'],transform=ax.transAxes)
    f.subplots_adjust(left=.03,right=.97,bottom=.03,top=.91)
    record_table('fig10_evidence_implications',[dict(evidence=r[0].replace('\n',' '),interpretation_limit=r[1].replace('\n',' '),development_implication=r[2].replace('\n',' ')) for r in content])
    assert status['status']['ROOT_CAUSE_PRIMARY']=='INCONCLUSIVE'
    return f

FUNCTIONS=[overview,activation,precursor,layer_forces,components,support,work_energy,positive_branch,stability,implications]
TRANSFORMS={
'fig01':'Read diagnostic onset, positive-state counts and normalized half-reactions. Boxes and arrows are an explanatory layout, not simulated data.',
'fig02':'strain*100; interval endpoints*100; active_count/bonds*100; median clamp depth=percentile index 2; lengths*1e6; depth/delta and depth/h.',
'fig03':'z_mid=(z_left+z_right)/2; z_mid/L; displacement_gap*1e6; engineering_strain*100. A_psi is the saved energy-density ratio; controls read directly from summary.',
'fig04':'z/L; saved layer_fz_over_R and cumulative_over_R. Top half only for cumulative external forces. R_cont=Pzz*L^2 per case.',
'fig05':'Bottom first-interior center particle; divide each saved signed algebraic contribution by its saved net axial force. Each pair sums to one.',
'fig06':'All archived particles (no subsampling, jitter or interpolation); force_norm_over_R versus d_over_delta, classified by full_support. Support counts and gates read directly.',
'fig07':'Absolute errors abs(energy_ratio-1) and abs(virtual_work_ratio-1); saved top_required ratio; archived 1% traction tolerance.',
'fig08':'Positive states maximum_s>1e-8; (strain-diagnostic_onset)*1e4 for phase/count; strain*100 for stress. Weighted mean uses existing pair quadrature weights, not phase-energy weighting.',
'fig09':'Stored physical Hessian/Schur minima; largest saved curvature_tol across FD steps per state. r_f/1e-9 and r_KKT/1e-8 from archived tolerances. Exact zero omitted, never clipped, on logarithmic axes.',
'fig10':'Actual onset relative shift, stable positive-state count, half-reaction ratios and support counts mapped to explicitly prospective implications.'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--only',nargs='*',help='Figure IDs for layout revision; omit to rebuild the complete collection')
    args=parser.parse_args();style.setup()
    for p in (FIG,OUT/'tables',MAN):p.mkdir(parents=True,exist_ok=True)
    write_plan(HERE)
    selected=[(s,fn) for s,fn in zip(SPECS,FUNCTIONS) if not args.only or s['id'] in args.only]
    # Freeze the small, declared source set before rendering. Never touch historical manifests.
    declared={p for s,_ in selected for p in s['sources']}
    before={r['path']:r for r in d.source_metadata(declared)}
    records=[];alltexts=[];allused=set(declared)
    pdf=PdfPages(OUT/'ALL_FIGURES.pdf') if not args.only else None
    for spec,fn in selected:
        d.USED.clear();fig=fn();used=set(spec['sources'])|d.USED;allused|=used
        # All implied sources must have been declared in the plan and frozen above.
        missing=used-declared
        if missing:raise RuntimeError(f'Undeclared source dependencies: {sorted(missing)}')
        stem=spec['id']+'_'+spec['slug'];fig.canvas.draw()
        texts=[x.get_text() for x in fig.findobj(Text) if x.get_text()]
        if any(re.search('[\u3400-\u9fff]',x) for x in texts):raise ValueError('Non-English figure text')
        alltexts.append('\n'.join([stem,*texts]))
        for ext in ('pdf','svg','png'):
            fig.savefig(FIG/(stem+'.'+ext),dpi=600,facecolor='white')
        if pdf:pdf.savefig(fig,facecolor='white')
        record=dict(id=spec['id'],title=spec['title'],files=[f'figures/{stem}.{ext}' for ext in ('pdf','svg','png')],
                    size_mm=(fig.get_size_inches()*25.4).tolist(),png_dpi=600,transformations=TRANSFORMS[spec['id']],
                    sources=d.source_metadata(used),caption=spec['caption'],no_new_simulation=True)
        (MAN/(spec['id']+'_sources.json')).write_text(json.dumps(record,indent=2)+'\n')
        records.append(record);plt.close(fig);print('Rendered',stem,flush=True)
    if pdf:pdf.close()
    after={r['path']:r for r in d.source_metadata(declared)}
    assert before==after,'A source file changed during rendering'
    if not args.only:
        (MAN/'figure_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
        (MAN/'selected_source_manifest.json').write_text(json.dumps(list(after.values()),indent=2)+'\n')
        (MAN/'rendered_figure_text.txt').write_text('\n\n'.join(alltexts)+'\n')
        captions=['# English figure captions','']
        for s in SPECS:
            stem=s['id']+'_'+s['slug'];captions.extend([f"## {s['id']} — {s['title']}",'',s['caption'],'',f'[PDF](output/figures/{stem}.pdf) · [SVG](output/figures/{stem}.svg) · [PNG](output/figures/{stem}.png)',''])
        (HERE/'CAPTIONS.md').write_text('\n'.join(captions))
        run=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__,
                 git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=d.ROOT,text=True).strip(),
                 figure_count=len(records),source_count=len(declared),sources_unchanged=True,
                 scientific_simulations_run=0,scientific_source_edits=0,scope='Archived-data plotting only; current worktree contains prior uncommitted changes.')
        (MAN/'RUN_VALIDATION.json').write_text(json.dumps(run,indent=2)+'\n')
        update_inventory(records)
    print('Complete: archived sources unchanged.',flush=True)

def update_inventory(records):
    path=MAN/'data_inventory.csv'
    if not path.exists():return
    with path.open(newline='') as f:rows=list(csv.DictReader(f));keys=list(rows[0])
    used={}
    for r in records:
        for src in r['sources']:used.setdefault(src['path'],[]).append(r['id'])
    for row in rows:
        if row['path'] in used:
            row['directly_usable']='yes_validated_for_listed_figures'
            row['inspection']='Loaded/validated and SHA256 recorded for '+', '.join(used[row['path']])
    with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

if __name__=='__main__':main()
