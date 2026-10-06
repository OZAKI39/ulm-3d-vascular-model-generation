"""Redraw exact saved result panels for readable projected labels; no new simulation."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
TASK='pfpd_baqp_finite_domain_boundary_consistency_v1'
V=ROOT/'verification'/TASK;T=ROOT/'results'/TASK/'tables'
SOURCES=[]
def rows(name):
 p=T/name;SOURCES.append(str(p.relative_to(ROOT)))
 with p.open() as f:return list(csv.DictReader(f))
def summary(c):
 p=V/c/'SUMMARY.json';SOURCES.append(str(p.relative_to(ROOT)));return json.loads(p.read_text())
def make():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'lines.linewidth':1.7,'lines.markersize':4.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','axes.axisbelow':True})
 out={};assets=HERE/'assets';assets.mkdir(exist_ok=True)
 def save(fig,name):
  for ext in ('png','pdf','svg'):fig.savefig(assets/f'{name}.{ext}',dpi=450,facecolor='white')
  out[name]=assets/f'{name}.png';plt.close(fig)
 fig,ax=plt.subplots(figsize=(4.05,2.10));fig.subplots_adjust(left=.19,right=.98,bottom=.27,top=.97)
 for c,color,marker in [('C5','#0072B2','o'),('C7','#D55E00','s')]:
  d=rows(f'{c}_layers.csv');ax.plot([float(r['z_m'])*1e6 for r in d],[float(r['layer_fz_over_R']) for r in d],marker=marker,color=color,label=c)
 ax.axhline(0,color='#777777',lw=.7);ax.set(xlabel='Reference z [μm]',ylabel='Layer force / Rcont',xticks=[-125,0,125],yticks=[-.8,0,.8],ylim=(-.88,.88));ax.legend(loc='upper right',frameon=False,ncol=2,handlelength=1.1,columnspacing=.7);ax.grid(axis='y',color='#e1e4e8',lw=.5);save(fig,'layer_forces')
 fig,ax=plt.subplots(figsize=(3.25,2.10));fig.subplots_adjust(left=.22,right=.97,bottom=.27,top=.97)
 d=[summary(c)['first_layer'][0] for c in ('C5','C7')];x=np.arange(2)
 ax.bar(x-.17,[r['direct_fraction'] for r in d],width=.34,color='#0072B2',label='Direct')
 ax.bar(x+.17,[r['Z_fraction'] for r in d],width=.34,color='#E69F00',label='Z correction')
 ax.set(xticks=x,xticklabels=['C5','C7'],ylabel='Fraction of net force',ylim=(0,1.2),yticks=[0,.5,1]);ax.legend(loc='upper center',ncol=2,frameon=False,handlelength=.8,columnspacing=.6);ax.grid(axis='y',color='#e1e4e8',lw=.5);save(fig,'direct_Z')
 fig,ax=plt.subplots(figsize=(3.75,2.10));fig.subplots_adjust(left=.23,right=.98,bottom=.27,top=.97)
 d=rows('L075_particles.csv');dep=np.array([float(r['d_over_delta']) for r in d]);val=np.array([float(r['force_norm_over_R']) for r in d]);full=np.array([r['full_support']=='True' for r in d]);assert full.sum()==summary('L075')['full_support_count']
 ax.scatter(dep[~full],val[~full],s=12,alpha=.4,color='#777777',marker='o',label='Truncated')
 ax.scatter(dep[full],val[full],s=25,alpha=.75,color='#D55E00',marker='D',label='Full family')
 ax.set_yscale('log');ax.axvline(1,ls='--',color='#555555',lw=.8)
 ax.set(xlabel='Nearest-surface depth / δ',ylabel='Nodal |f| / Rcont',xticks=[0,.5,1,1.5],yticks=[1e-16,1e-10,1e-4],ylim=(1e-18,1e1),xlim=(-.06,1.52));ax.legend(loc='upper center',frameon=False,ncol=2,handletextpad=.2,columnspacing=.4,markerscale=.9);ax.grid(axis='y',color='#e1e4e8',lw=.5);save(fig,'full_family')
 (assets/'inset_sources.json').write_text(json.dumps({'source_files':sorted(set(SOURCES)),'method':'Exact existing layer CSV, nodal CSV and algebraic summary values; no interpolation, filtering of observations, solver call or fabricated trend. The same panels as the existing audit figures are re-typeset at readable inset font sizes.','source_figure_counterparts':[f'results/{TASK}/figures/'+n+'.png' for n in ['affine_internal_force_layers','first_layer_force_decomposition','large_domain_affine_patch_test']]},indent=2)+'\n')
 return out
if __name__=='__main__':make()
