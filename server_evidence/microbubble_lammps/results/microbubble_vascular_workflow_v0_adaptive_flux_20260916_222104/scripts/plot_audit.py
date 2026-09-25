from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import Patch,FancyBboxPatch
from flow_geometry import FrozenSampler
S=Path(__file__).resolve().parents[1];sec=np.load(S/'geometry/INJECTION_SECTION.npz');a=json.loads((S/'validation/INLET_FLUX_AUDIT.json').read_text());f=FrozenSampler(S/'fields/FROZEN_FLOW_FIELD_V0.h5')
axes=np.c_[sec['axis1'],sec['axis2']];xy=(sec['points']-sec['center'])@axes;bounds=np.c_[xy.min(axis=0),xy.max(axis=0)]
x=np.linspace(*bounds[0],301);y=np.linspace(*bounds[1],301);xx,yy=np.meshgrid(x,y);coords=np.c_[xx.ravel(),yy.ravel()]
inside=MplPath(xy).contains_points(coords);p=sec['center']+coords@axes.T;st,u=f.query(p);valid=inside&(st==0);normal_speed=u@(-sec['normal'])
plt.rcParams.update({'font.size':10,'savefig.dpi':170})
def save(name):plt.savefig(S/'visualization'/name,bbox_inches='tight');plt.close()
fig,ax=plt.subplots(figsize=(7,6));C=np.full(len(coords),np.nan);C[valid]=np.maximum(normal_speed[valid],0)*1e3
im=ax.pcolormesh(xx*1e6,yy*1e6,C.reshape(xx.shape),shading='auto',cmap='viridis');plt.colorbar(im,ax=ax,label='Known positive normal velocity (mm/s)')
invalid=(inside&~valid).reshape(xx.shape);ax.contourf(xx*1e6,yy*1e6,invalid.astype(float),levels=[.5,1.5],colors=['#df7290'],alpha=.9)
loop=np.vstack([xy,xy[0]])*1e6;ax.plot(*loop.T,color='black',lw=1);ax.set_aspect('equal');ax.set_xlabel('Plane tangent 1 (µm)');ax.set_ylabel('Plane tangent 2 (µm)');ax.legend(handles=[Patch(facecolor='#df7290',label='Undefined: native sampler SOLID status')],loc='upper left',bbox_to_anchor=(0,-.13));ax.set_title('Real inlet velocity coverage: incomplete\nInvalid values are not zero velocity');save('inlet_positive_flux_map.png')
# Alternate requested map name contains the same diagnostic, not a second measurement.
(S/'visualization/inlet_flux_map.png').write_bytes((S/'visualization/inlet_positive_flux_map.png').read_bytes())
r=a['quadrature'];n=[x['quadrature_points'] for x in r];qp=[x['partial_Q_positive_m3_s']*1e15 for x in r];coverage=[x['valid_area_fraction']*100 for x in r]
fig,axs=plt.subplots(1,2,figsize=(11,4.5));axs[0].plot(n,qp,'o-');axs[0].set_xscale('log');axs[0].set_xlabel('Triangle quadrature points');axs[0].set_ylabel('Partial positive flux (pL/s)');axs[0].set_title('Defined portion only; full Q remains unknown')
axs[1].plot(n,coverage,'o-');axs[1].axhline(100,color='grey',ls='--');axs[1].set_xscale('log');axs[1].set_ylim(75,101);axs[1].set_xlabel('Triangle quadrature points');axs[1].set_ylabel('Area with valid velocity (%)');axs[1].set_title('Refinement does not fill missing support');fig.suptitle('Partial integral convergence does not establish full inlet flux');fig.tight_layout();save('inlet_flux_convergence.png')
fig,ax=plt.subplots(figsize=(10,10));ax.set_xlim(0,10);ax.set_ylim(0,12);ax.axis('off')
nodes=[('Frozen field + real inlet section','inherited, hash verified'),('Full positive inlet volume flux Q(t)','BLOCKED: 16.4% of area has undefined velocity'),('Concentration or target number flux','production settings UNSPECIFIED'),('Flux accumulator → source SonoVue draw','planned; not implemented after STOP'),('Size admissibility → position search','planned; ADMIT / PENDING / SIZE REJECT'),('LAMMPS + existing rigid-sphere RK2','planned; no workflow trajectory run'),('Open-port exits + particle / flux accounting','planned; no lifecycle validation run')]
for i,(title,detail) in enumerate(nodes):
    yy=10.6-i*1.5;color='#fbe0e5' if i==1 else ('#dcecdc' if i==0 else '#f0f0f0');edge='#b2182b' if i==1 else '#666666'
    ax.add_patch(FancyBboxPatch((.45,yy),9.1,1.05,boxstyle='round,pad=.12',facecolor=color,edgecolor=edge,lw=1.5))
    ax.text(5,yy+.70,title,ha='center',va='center',fontsize=11);ax.text(5,yy+.28,detail,ha='center',va='center',fontsize=9,color=edge)
    if i<len(nodes)-1:ax.annotate('',xy=(5,yy-.36),xytext=(5,yy-.13),arrowprops={'arrowstyle':'->','color':'grey','linestyle':'--' if i>=1 else '-'})
ax.set_title('Adaptive injection workflow: implementation stopped at flux audit',pad=20);save('adaptive_injection_workflow.png')
print('Diagnostic flux/coverage/convergence and explicitly planned workflow figures written; no particle data fabricated.')
