"""Plot read-only internal flux checks, preserving the solved velocity field."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/internal_sections'
rows=list(csv.DictReader((OUT/'section_flux.csv').open()))
div=json.loads((OUT/'DIVERGENCE_THEOREM_CHECK.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
colors=['#355c9c','#bd5e18','#319678','#9866b0']
roles=['INLET','OUTLET_01','OUTLET_02','OUTLET_03']
for role,color in zip(roles,colors):
    subset=[r for r in rows if r['boundary']==role]
    x=[0]+[float(r['inward_offset_mm']) for r in subset]
    y=[100]+[100*abs(float(r['signed_outward_Q_m3_s'])/float(r['signed_boundary_reference_m3_s'])) for r in subset]
    axes[0].plot(x,y,'o-',color=color,label=role)
axes[0].axhline(100,color='gray',ls='--',lw=.7)
axes[0].set(xlabel='Distance inward from port (mm)',ylabel='Internal section flow / end-face flow (%)',title='Same frozen P1 velocity; no flow correction')
axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
x=np.arange(4)
q=[r['surface_flux_difference_m3_s']*6e10 for r in div['rows']]
d=[r['volume_divergence_integral_m3_s']*6e10 for r in div['rows']]
axes[1].bar(x-.17,q,width=.34,label='End face minus 2-mm section')
axes[1].bar(x+.17,d,width=.34,label='Independent volume integral of div(u)')
axes[1].set(xticks=x,xticklabels=roles,ylabel='Signed flux difference (μL/min)',title='Divergence theorem cross-check of 2-mm slabs')
axes[1].tick_params(axis='x',labelsize=8);axes[1].legend(fontsize=8);axes[1].grid(axis='y',alpha=.2)
fig.savefig(OUT/'internal_flux_and_divergence.png',dpi=300)
fig.savefig(OUT/'internal_flux_and_divergence.pdf')
plt.close(fig)
print('Saved',OUT/'internal_flux_and_divergence.png')
