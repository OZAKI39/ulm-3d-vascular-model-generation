"""Static review figures and their CSV sources; no visual-only acceptance gates."""
import json
import os
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR','/tmp/cf2003_mplconfig_20260916_192326')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120})
data=np.genfromtxt(ROOT/'raw/CF2003_REFERENCE_DENSE_GRID.csv',delimiter=',',names=True)
points=np.genfromtxt(ROOT/'reference/CF2003_FREE_SHEAR_REFERENCE_V0.csv',delimiter=',',names=True)
x=data['epsilon']
def save(fig,name):
    fig.tight_layout();fig.savefig(ROOT/'visualization'/name,dpi=170);plt.close(fig)
def axis(ax,title,ylabel):
    ax.set(xscale='log',xlim=(.001,.2),xlabel='Surface gap / radius, epsilon = h/a',ylabel=ylabel,title=title)
    ax.grid(alpha=.2)
for key,name,ylabel,color in [('FU','FU_CF2003_vs_epsilon.png','FU = Ux / (kappa*l)','#184b87'),
                             ('FOMEGA','FOMEGA_CF2003_vs_epsilon.png','FOMEGA = Omega_y / (kappa/2)','#117b60')]:
    fig,ax=plt.subplots(figsize=(7.4,4.5));ax.plot(x,data[key],color=color,label='CF2003 Table 17; checked 2012 correction')
    ax.scatter(points['epsilon'],points[key],s=20,color=color,zorder=3,label='Required reference grid')
    axis(ax,'Certified source-fit domain: 0.001 <= epsilon <= 0.2',ylabel);ax.legend(fontsize=8);save(fig,name)
fig,ax=plt.subplots(figsize=(7.4,4.5));ax.plot(x,data['Omega_a_over_U'],color='#6d438c',label='CF2003: FOMEGA / [2(1+epsilon) FU]')
ax.axhline(.5676,color='gray',ls='--',label='Historical contact limit (epsilon -> 0), outside domain')
axis(ax,'Rolling ratio: finite gap does not equal contact','Omega*a / U (positive right-hand convention)');ax.legend(fontsize=8);save(fig,'rolling_ratio_vs_epsilon.png')
fig,axes=plt.subplots(1,2,figsize=(11.2,4.5))
for ax,key,old,label in zip(axes,['FU','FOMEGA'],['GCB_old_FU','GCB_old_FOMEGA'],['Translation','Rotation']):
    ax.plot(x,data[key],lw=2,color='#184b87',label='CF2003 primary candidate (source validated)')
    ax.plot(x,data[old],ls='--',color='#b36128',label='GCB old expression: historical diagnostic')
    if key=='FU':ax.plot(x,data['GCB_corrected_FU'],ls=':',color='#587957',label='Corrected first order: cross-check only')
    ax.scatter([.001],[.3966 if key=='FU' else .4268],marker='x',color='black',label='GCB Table 3 transcription')
    axis(ax,label,key);ax.legend(fontsize=7)
save(fig,'GCB_vs_CF2003_diagnostic.png')
table=Path('/home/lzy/projects/compre_output/wall_hydrodynamics_v0/20260916_130547/tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5')
with h5py.File(table,'r') as f:
    eps=f['epsilon'][:];r=f['R_total_scaled'][:];sel=eps<=.2;eps=eps[sel];r=r[sel]
np.savetxt(ROOT/'raw/RESISTANCE_COMPONENTS.csv',np.c_[eps,r[:,2,2],r[:,0,0],r[:,4,4],r[:,0,4]],delimiter=',',header='epsilon,normal_TT,tangential_TT,parallel_RR,TR_scaled',comments='')
fig,axes=plt.subplots(1,2,figsize=(10.5,4.4))
axes[0].loglog(eps,r[:,2,2],label='Normal TT / bulk TT');axes[0].loglog(eps,r[:,0,0],label='Tangential TT / bulk TT')
axis(axes[0],'Frozen RMBW total resistance','Dimensionless resistance');axes[0].legend(fontsize=8)
axes[1].semilogx(eps,r[:,4,4],label='Parallel RR / bulk RR');axes[1].semilogx(eps,r[:,0,4],label='TR / sqrt(bulk TT * bulk RR)')
axis(axes[1],'RR/TR: inherited finite-gap accuracy limitations','Dimensionless resistance');axes[1].legend(fontsize=8)
save(fig,'resistance_vs_epsilon.png')
shear=np.genfromtxt(ROOT/'raw/SHEAR_RHS_COMPONENTS.csv',delimiter=',',names=True)
fig,axes=plt.subplots(1,3,figsize=(12.4,4.1))
for key in ['Fx_N','Fy_N','Fz_N']:axes[0].semilogx(shear['epsilon'],shear[key],label=key)
for key in ['Tx_Nm','Ty_Nm','Tz_Nm']:axes[1].semilogx(shear['epsilon'],shear[key],label=key)
axes[2].semilogx(shear['epsilon'],shear['Fz_N'],color='#a23357',label='Normal force Fz')
axis(axes[0],'Analytic shear RHS force','N');axis(axes[1],'Analytic shear RHS torque','N m')
axis(axes[2],'Normal force; max abs = %.3g N'%np.max(abs(shear['Fz_N'])),'N')
axes[2].set_ylim(-1e-24,1e-24)
for ax in axes:ax.legend(fontsize=8)
fig.suptitle('a=1 micrometer, mu=0.001 Pa s, g=(100,0,0) 1/s, n=+z',fontsize=11)
fig.tight_layout(rect=(0,0,1,.94));fig.savefig(ROOT/'visualization/shear_rhs_components_vs_epsilon.png',dpi=170);plt.close(fig)
print('Saved 6 scientific figures; corresponding CSV data retained')
