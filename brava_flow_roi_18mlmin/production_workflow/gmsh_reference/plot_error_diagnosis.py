"""Exactly two diagnostic figures, each PNG 300 dpi and vector PDF."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'gmsh_reference/error_diagnosis'
FIG=ROOT/'gmsh_reference/figures'
def read(n):return json.loads((ROOT/n).read_text())
def main():
 d=read('V2_FIELD_ERROR_DECOMPOSITION.json');p=read('P1_FIELD_PARTICLE_SENSITIVITY.json')
 plt.rcParams.update({'font.size':12})
 fig,axes=plt.subplots(2,2,figsize=(14,10),constrained_layout=True)
 for ax,role,title in zip(axes[0],['PARTICLE_RELEVANT','BULK'],['Particle velocity bands','Bulk velocity bands']):
  rows=[r for r in d['groups'] if r['role']==role];x=np.arange(len(rows))
  for offset,key,label,color in [(-.24,'space_relative_L2','P1 space','#b7812b'),(0,'discrete_relative_L2','CFD minus P1','#5788a5'),(.24,'total_relative_L2','CFD total','#bf4141')]:
   ax.bar(x+offset,[r['velocity'][key]*100 for r in rows],width=.24,label=label,color=color)
  ax.axhline(5,color='black',ls='--',lw=1,label='Frozen 5% gate')
  ax.set_xticks(x,[f"{r['wall_distance_m']*1e6:.3g}" for r in rows]);ax.set(title=title,xlabel='Wall distance (µm)',ylabel='Relative L2 (%)')
  ax.legend(fontsize=9,ncol=2);ax.grid(axis='y',alpha=.2)
 ax=axes[1,0];profile=d['radial_profile'];depth=np.array(profile['wall_distances_m'])*1e6
 for key,label,color in [('P1','P1 analytic interpolant','#b7812b'),('CFD','Frozen CFD','#bf4141')]:
  ax.semilogx(depth,np.array(profile['fields'][key]['relative_axial_L2'])*100,label=label,color=color)
 for value in [3,4,5]:ax.axhline(value,ls='--',color='gray',lw=.8);ax.text(depth[-1]*1.1,value,str(value)+'%',fontsize=9)
 ymax=max(max(v['relative_axial_L2']) for v in profile['fields'].values())*105
 ax.set(title='Central radial profile: 16 fixed rays',xlabel='Wall distance (µm)',ylabel='Axial relative L2 (%)',ylim=(0,ymax),xlim=(.1,900));ax.legend(fontsize=10);ax.grid(alpha=.2)
 ax=axes[1,1];rows=[r for r in d['groups'] if r['role']!='WALL_STRESS_DIAGNOSTIC'];depth=np.array([r['wall_distance_m'] for r in rows])*1e6
 for name,key,color in [('CFD total','total_absolute','#bf4141'),('P1 space','space_absolute','#b7812b'),('CFD minus P1','discrete_absolute','#5788a5')]:
  ax.loglog(depth,[r['velocity'][key]['RMS_mm_s'] for r in rows],'o-',label=name,color=color)
 ax.set(title='Absolute velocity error at frozen bands',xlabel='Wall distance (µm)',ylabel='RMS error (mm/s)');ax.legend(fontsize=10);ax.grid(alpha=.2)
 fig.suptitle('Frozen V2 velocity error decomposition — diagnostic only\nVector norms do not add; all cross terms and original gates are retained',fontsize=16)
 for ext in ['png']:fig.savefig(FIG/f'v2_error_decomposition.{ext}',dpi=300)
 plt.close(fig)
 fig,axes=plt.subplots(3,4,figsize=(19,11),constrained_layout=True)
 specs=[('A',.375,.002),('B',.964365,.002),('C',.964365,2),('D',.964365,10)]
 for col,(name,radius,gap) in enumerate(specs):
  with np.load(OUT/f'comparison_{name}.npz') as z:data={k:z[k] for k in z.files}
  t=data['time_s']*1e3;metric=next(m for m in p['comparison'] if m['case']==name)
  ax=axes[0,col]
  ax.plot(t,data['exact_speed_m_s']*1e3,label='Exact',color='#246ba0');ax.plot(t,data['P1_speed_m_s']*1e3,label='P1 analytic',color='#b7812b')
  ax.set(title=f'{name}: a={radius:g} µm, gap={gap:g} µm',ylabel='Particle RHS speed (mm/s)',xlabel='Physical time (ms)');ax.legend(fontsize=10,loc='upper right',bbox_to_anchor=(.98,.65))
  ax.text(.04,.08,f"Speed L2: {metric['speed_history_relative_L2']*100:.3f}%\nTransit error: {metric['transit_relative_error']*100:.3f}%",transform=ax.transAxes,fontsize=10)
  ax=axes[1,col]
  for kind,label,color in [('exact','Exact','#246ba0'),('P1','P1 analytic','#b7812b')]:
   pos=data[kind+'_position_m'];ax.plot(t,(pos[:,2]-pos[0,2])*1e6,label=label,color=color)
  ax.set(xlabel='Physical time (ms)',ylabel='Axial displacement (µm)')
  ax=axes[2,col];exact_speed=data['exact_speed_m_s'];diff=(data['P1_speed_m_s']/exact_speed-1)*100
  ax.plot(t,abs(diff),color='#bf4141',label='Absolute speed difference')
  ax.axhline(2,color='gray',ls='--',label='2% reference');ax.axhline(5,color='black',ls=':',label='5% reference')
  ax.set(xlabel='Physical time (ms)',ylabel='Speed difference (%)',ylim=(0,6.5))
  ax.text(.04,.98,'Radial separation: 0\nGap difference: 0\nContact category unchanged',transform=ax.transAxes,va='top',fontsize=10)
  if col==0:ax.legend(loc='lower right',fontsize=8)
 for ax in axes.flat:ax.grid(alpha=.2)
 fig.suptitle('Exact vs same-mesh P1 analytic field — production particle dynamics\nDiagnostic only: dt max 0.125 ms; 0.2 s; 10 µm local transit; no CFD trajectories',fontsize=16)
 for ext in ['png']:fig.savefig(FIG/f'p1_particle_sensitivity.{ext}',dpi=300)
 plt.close(fig)
if __name__=='__main__':main()
