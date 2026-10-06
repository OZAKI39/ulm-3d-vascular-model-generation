"""Vector-friendly 180 mm paper theme; Okabe-Ito palette and redundant markers."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
MM=1/25.4
COLORS={'C5':'#0072B2','C7':'#D55E00','L050':'#009E73','L075':'#CC79A7','direct':'#0072B2','Z':'#E69F00','gray':'#777777','ink':'#222222'}
MARKERS={'C5':'o','C7':'s','L050':'^','L075':'D'}
def setup():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'axes.titleweight':'normal','xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'figure.titlesize':10,'axes.linewidth':.65,'lines.linewidth':1.15,'lines.markersize':4,'patch.linewidth':.65,'xtick.major.width':.6,'ytick.major.width':.6,'xtick.direction':'out','ytick.direction':'out','axes.spines.top':False,'axes.spines.right':False,'axes.axisbelow':True,'grid.color':'#dddddd','grid.linewidth':.45,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':600,'figure.dpi':110,'axes.formatter.use_mathtext':True})
def panels(rows=1,cols=2,height=85,width=180):
 f,axs=plt.subplots(rows,cols,figsize=(width*MM,height*MM),squeeze=False,layout='constrained')
 f.set_constrained_layout_pads(w_pad=.06,h_pad=.05,wspace=.08,hspace=.08)
 for i,ax in enumerate(axs.flat):
  ax.text(-.12,1.075,f'({chr(97+i)})',transform=ax.transAxes,fontweight='bold',fontsize=9,va='bottom')
  ax.grid(axis='y',alpha=.65)
 return f,axs.flatten()
def legend(ax,**kw):return ax.legend(frameon=False,**kw)
