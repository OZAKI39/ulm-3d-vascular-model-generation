"""Actual 2.0 ms components, shown before interpreting the full protocol."""
from common import *
from figures import plt,ax3,Line3DCollection,COLORS,NAMES,save
from multifragment_shear_objective import shape

def main():
 ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'BIPHASIC_FIELD.npz');graph=module('rise_end_figure_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py');records=[];width=0.
 for mode in CASES:
  p=O/'cases'/mode/'keys/rise_end.npz';s=load_npz(p);assert abs(float(s['time'])-.002)<1e-14
  lab,sz,att=graph.topology(ref,s['d']);free_labels=np.flatnonzero(~att);label=free_labels[np.argmax(sz[free_labels])];ids=np.flatnonzero(lab==label)
  live=np.isin(ref['bonds'],ids).all(1)&(s['d']<1);outgoing=(np.isin(ref['bonds'],ids).sum(1)==1)&(s['d']<1);assert not outgoing.any()
  record=dict(case=mode,state=s,source=str(p.relative_to(O)),source_sha256=sha(p),ids=ids,bonds=ref['bonds'][live],reference_shape=shape(ref['X'][ids]),current_shape=shape(s['x'][ids]),strong_reference_nodes=int((field['material_phase'][ids]==1).sum()))
  records.append(record);width=max(width,1.25*np.ptp(ref['X'][ids]*1e6,axis=0).max(),1.25*np.ptp(s['x'][ids]*1e6,axis=0).max())
 fig=plt.figure(figsize=(15,12));fig.subplots_adjust(wspace=.18,hspace=.34,bottom=.11,top=.85)
 for row,r in enumerate(records):
  ids=r['ids'];mode=r['case'];s=r['state']
  for col,(xyz,label) in enumerate([(ref['X'],'Reference positions of the same original IDs'),(s['x'],'Actual accepted positions at 2.000 ms')]):
   center=(xyz[ids].min(0)+xyz[ids].max(0))*.5e6
   ax=ax3(fig,(2,2,row*2+col+1),NAMES[mode]+f' · {len(ids)} nodes\n'+label,(center-width/2,center+width/2))
   dense=len(r['bonds'])>1000
   ax.add_collection3d(Line3DCollection(xyz[r['bonds']]*1e6,colors='#244e63',linewidths=.3 if dense else .7,alpha=.12 if dense else .55,zorder=1))
   ax.scatter(*xyz[ids].T*1e6,s=22 if dense else 28,c=COLORS[mode],depthshade=False,zorder=3)
   ax.text2D(.02,.96,f"{len(r['bonds'])} surviving internal bonds; zero live outgoing bonds\n"+('20–200 count band; elongated reference shape' if len(ids)<=200 else '>200 band; not counted as a target fragment'),transform=ax.transAxes,fontsize=10,va='top')
 fig.suptitle('Rise-end stage comparison · actual original-bond graph components\nReference and deformed positions · same camera and local physical width · no artificial separation',fontsize=17)
 fig.text(.05,.035,f'Displayed width: {width+2*max(.06*width,20.):.1f} µm in every panel. Absolute physical coordinates; displacement ×1.\nSurviving edges are from the 2.0 ms accepted state in both reference and deformed views. This stage is not the final outcome.',fontsize=12)
 save(fig,'FIG_07E_RISE_END_REFERENCE_CURRENT')
 evidence=[]
 for r in records:
  evidence.append({k:v for k,v in r.items() if k not in ['state','ids','bonds']}|dict(original_atom_ids=r['ids']+1,internal_live_bonds=len(r['bonds']),time_s=float(r['state']['time'])))
 dump(O/'RISE_END_COMPONENT_FIGURE_EVIDENCE.json',dict(records=evidence,display_width_um=width+2*max(.06*width,20.),scope='Same 2.0 ms accepted stage, not final protocol outcome. U strong_reference_nodes are spatial covariates only.'))
 print('RISE_END_COMPONENT_FIGURE_DONE',[(r['case'],len(r['ids']),len(r['bonds'])) for r in records],flush=True)

if __name__=='__main__':main()
