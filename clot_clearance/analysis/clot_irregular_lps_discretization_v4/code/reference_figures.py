"""Measured geometry and elastic figures, independent of the formal outcome."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import shutil

def main():
 out=F/'figures';out.mkdir(parents=True,exist_ok=True);old=load_npz(V3_O/'provenance/reference.npz');new=load_npz(O/'provenance/reference.npz');field=load_npz(O/'BIPHASIC_FIELD.npz');n=int(new['clot_count'])
 plt.rcParams.update({'font.size':13,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#a8b0b7','grid.alpha':.65,'savefig.facecolor':'white'})
 fig=plt.figure(figsize=(15,10),layout='constrained')
 subset=(abs(old['X'][:,0])<.18e-3)&(abs(old['X'][:,1])<.18e-3)&(old['X'][:,2]>.075e-3)&(old['X'][:,2]<.3e-3)
 for col,(name,ref) in enumerate([('V3 regular',old),('V4 irregular',new)]):
  for row in [0,1]:
   ax=fig.add_subplot(2,2,1+col+2*row,projection='3d');ids=np.arange(len(ref['X'])) if row==0 else np.flatnonzero(subset);p=ref['X'][ids]*1e3;colors=np.where(ids>=n,'#7499ad',np.where(field['material_phase'][ids]==1,'#6f1025','#c64b5a'))
   ax.scatter(*p.T,s=3 if row==0 else 26,c=colors,alpha=.65 if row==0 else .95,depthshade=True,linewidths=0)
   if row:
    b=ref['bonds'];mask=subset[b].all(1)&~ref['initial_cut'];from mpl_toolkits.mplot3d.art3d import Line3DCollection
    ax.add_collection3d(Line3DCollection(ref['X'][b[mask]]*1e3,colors='#aeb8c0',linewidths=.2,alpha=.2))
   ax.set_title(name+(' — full reference' if row==0 else ' — fixed local window'))
   ax.set_xlabel('X (mm)',labelpad=10);ax.set_ylabel('Y (mm)',labelpad=10);ax.set_zlabel('Z (mm)',labelpad=10);ax.view_init(elev=15,azim=-58,roll=0)
   bounds=[[-.65,.65],[-.65,.65],[-.5,.35]] if row==0 else [[-.18,.18],[-.18,.18],[.075,.3]]
   ax.set_xlim(bounds[0]);ax.set_ylim(bounds[1]);ax.set_zlim(bounds[2]);ax.set_box_aspect(np.ptp(bounds,axis=1))
 fig.suptitle('Same material IDs, boundary, phases and nominal scale\nReference positions, original bonds and quadrature weights change',fontsize=17)
 fig.savefig(out/'FIG_03_REFERENCE_POINT_CLOUDS.png',dpi=200);plt.close(fig);shutil.copy2(out/'FIG_03_REFERENCE_POINT_CLOUDS.png',out/'FIG_03_REGULAR_VS_IRREGULAR_LATTICE.png')
 fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained');stats={}
 for name,ref,color in [('V3 regular',old,'#718c9e'),('V4 irregular',new,'#aa203e')]:
  b=ref['bonds'];bulk=ref['bond_types']==1;dirs=(ref['X'][b[bulk,1]]-ref['X'][b[bulk,0]])/ref['r0'][bulk,None];neighbors=np.bincount(b.ravel(),minlength=len(ref['X']))[:n]
  axes[0,0].hist(neighbors,bins=np.arange(20,131,3),histtype='step',density=True,color=color,label=name,lw=2)
  axes[0,1].hist(ref['r0'][bulk]*1e6,bins=np.linspace(25,151,55),histtype='step',density=True,color=color,lw=2)
  axes[1,0].hist(abs(dirs).ravel(),bins=np.linspace(0,1,41),histtype='step',density=True,color=color,lw=2)
  axes[1,1].hist(ref['V'][:n]/V,bins=np.linspace(.8,1.2,45),histtype='step',density=True,color=color,lw=2)
  stats[name]=dict(bonds=len(b),material_phase_bond_counts={str(k):int(np.sum(bulk&(field['material_phase'][b[:,0]]+field['material_phase'][b[:,1]]==k))) for k in [0,1,2]},direction_tensor=dirs.T@dirs/len(dirs),axis_alignment_fourth_moment=float(np.mean(np.sum(dirs**4,axis=1))),neighbor_min_mean_max=[int(neighbors.min()),float(neighbors.mean()),int(neighbors.max())])
 for ax,title,xlabel in zip(axes.ravel(),['Reference neighbor counts','Bulk reference bond lengths','Direction components','Material quadrature volumes'],['Neighbors / material node','Length (µm)','|n·coordinate axis|','Vi / (50 µm)³']):ax.set_title(title);ax.set_xlabel(xlabel);ax.set_ylabel('Probability density')
 axes[0,0].legend();fig.savefig(out/'FIG_04_NEIGHBORHOOD_AND_VOLUME.png',dpi=190);plt.close(fig);dump(O/'REFERENCE_DIRECTION_STATISTICS.json',stats)
 q=json.loads((O/'ELASTIC_PATCH_CHECK.json').read_text());fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
 for j,name in enumerate(['regular','irregular']):
  values=[r['modulus_energy_Pa'] for r in q['patches'][name]['rows']];color=['#718c9e','#aa203e'][j]
  axes[0].bar(j,values[0]/1e3,color=color);axes[1].bar(np.arange(3)+(j-.5)*.35,values[1:],width=.35,color=color,label=name)
 axes[0].set_xticks([0,1],['Regular','Irregular']);axes[0].set_ylabel('Affine bulk response (kPa)');axes[0].set_title('Volume strain ε = 10⁻⁴')
 axes[1].set_xticks(np.arange(3),['XY','XZ','YZ']);axes[1].set_ylabel('Affine shear response (Pa)');axes[1].set_title('Simple shear γ = 10⁻⁴');axes[1].legend()
 changes=q['relative_modulus_changes'];axes[2].bar(list(changes),100*np.array(list(changes.values())),color='#aa203e');axes[2].axhline(15,color='#333',ls='--');axes[2].axhline(-15,color='#333',ls='--');axes[2].set_ylim(-17,17);axes[2].set_ylabel('V4 / V3 − 1 (%)');axes[2].set_title('Frozen ±15% comparison gate')
 fig.suptitle('729-node affine patches; prescribed coordinates, no material refit',fontsize=16);fig.savefig(out/'FIG_05_ELASTIC_RESPONSE.png',dpi=190);plt.close(fig)
 # Equal-volume limit of the generalized full endpoint-star energy, including
 # partially and completely damaged bonds, without another dynamics case.
 r=load_npz(O/'tests/elastic/regular/provenance/reference.npz');m0=module('old_energy',R/'analysis/softening_feasibility_v1/model.py').SofteningLPS(r);m1=energy_model(r);rng=np.random.default_rng(711);d=rng.random(len(r['bonds']));d[::17]=1;d[::11]=0;x=r['X']@np.array([[1.02,.01,0],[0,1.01,.005],[.003,0,.99]]).T;z0=m0.evaluate(x,d);z1=m1.evaluate(x,d)
 errors=dict(force_max_N=float(abs(z0['force']-z1['force']).max()),energy_J=abs(z0['U']-z1['U']));assert errors['force_max_N']<1e-14 and errors['energy_J']<1e-20;dump(O/'UNIFORM_VOLUME_LIMIT_CHECK.json',dict(status='PASS',errors=errors,random_damage_is_offline_algebra_fixture=True))
 print('REFERENCE_FIGURES_PASS',errors,flush=True)
if __name__=='__main__':main()
