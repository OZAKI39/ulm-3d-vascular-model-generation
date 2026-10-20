"""Read-only 0.8 ms snapshot comparison; never a final protocol conclusion."""
from common import *
from figures import plt,ax3,scatter,Line3DCollection,DAMAGE,ScalarMappable,Normalize

def main():
 ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'BIPHASIC_FIELD.npz');b=ref['bonds'];states={};evidence={};target=.0008
 for mode in CASES:
  # Regular snapshots are atomically published and never overwritten.
  candidates=[]
  for p in sorted((O/'cases'/mode/'keys').glob('frame_*.npz')):
   if '.tmp.' in p.name:continue
   with np.load(p,allow_pickle=False) as z:t=float(z['time'])
   candidates.append((abs(t-target),p,t))
  _,p,t=min(candidates);assert abs(t-target)<1e-6
  s=load_npz(p);states[mode]=s;bulk=(ref['bond_types']==1)&~ref['initial_cut'];groups={}
  for name,sel in [('WEAK_WEAK',bulk&(field['material_phase'][b]==0).all(1)),('WEAK_STRONG',field['cross_phase']&~ref['initial_cut']),('STRONG_STRONG',bulk&(field['material_phase'][b]==1).all(1)),('INTERFACE',ref['bond_types']==2)]:
   groups[name]=dict(eligible=int(sel.sum()),softened=int((sel&(s['d']>0)).sum()),completely_failed=int((sel&(s['d']==1)).sum()))
  evidence[mode]=dict(source=str(p.relative_to(O)),sha256=sha(p),actual_time_s=t,phase_pair_counts=groups,scope='U labels are H reference spatial covariates only; U uses uniform thresholds everywhere.')
 limits=(np.concatenate([s['x'] for s in states.values()]).min(0)*1e6,np.concatenate([s['x'] for s in states.values()]).max(0)*1e6)
 fig=plt.figure(figsize=(15,8));fig.subplots_adjust(wspace=.18,bottom=.2,top=.83)
 for j,mode in enumerate(CASES):
  s=states[mode];ax=ax3(fig,(1,2,j+1),('U · Uniform' if mode=='uniform' else 'H · Biphasic')+f' | actual t = {float(s["time"])*1e3:.5f} ms',limits);scatter(ax,ref,s,size=8)
  sel=(s['d']==1)&~ref['initial_cut']&(ref['r0']<=1.01*DX)
  if sel.any():ax.add_collection3d(Line3DCollection(s['x'][b[sel]]*1e6,colors='#9f254c',linewidths=.65,alpha=.55))
  g=evidence[mode]['phase_pair_counts'];ax.text2D(.01,.01,'Natural full failures: '+str(sum(v['completely_failed'] for v in g.values()))+'\nFull protocol comparison pending',transform=ax.transAxes,fontsize=12)
 cax=fig.add_axes([.31,.1,.38,.025]);fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap=DAMAGE),cax=cax,orientation='horizontal',label='Mean incident damage d (initial notch excluded)')
 fig.suptitle('Early accepted-state comparison · same camera and physical scale · displacement ×1\nDark traces: naturally failed nearest bonds; phase labels do not define fragments',fontsize=17)
 out=F/'previews';out.mkdir(exist_ok=True);fig.savefig(out/'UH_EARLY_DAMAGE_0P8MS.png',dpi=240,bbox_inches='tight');fig.savefig(out/'UH_EARLY_DAMAGE_0P8MS.pdf',bbox_inches='tight');plt.close(fig)
 dump(O/'UH_EARLY_DAMAGE_0P8MS_READBACK.json',dict(status='READ_ONLY_INTERIM_COMPARISON',cases=evidence));print(json.dumps(evidence,default=cv),flush=True)

if __name__=='__main__':main()
