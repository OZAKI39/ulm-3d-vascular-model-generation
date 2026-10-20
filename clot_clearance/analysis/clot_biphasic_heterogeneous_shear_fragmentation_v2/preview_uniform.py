"""Completed uniform-control readback only; no prediction of the biphasic run."""
from figures import *

ref=load_npz(O/'provenance/reference.npz');case=O/'cases/uniform';summary=json.loads((case/'SUMMARY.json').read_text())
assert summary['status']=='COMPLETED'
s=load_npz(case/'keys/final.npz');lab,sz,att=graph.topology(ref,s['d']);b=ref['bonds']
registry=json.loads((case/'COMPONENT_REGISTRY.json').read_text());lookup={tuple(q['atom_ids']):int(cid) for cid,q in registry.items()}
labels=np.flatnonzero(~att&(sz>=20)&(sz<=200));labels=labels[np.argsort(-sz[labels])];records=[]
assert len(labels)>0,'No target component exists to draw'
ids_all=[np.flatnonzero(lab==label) for label in labels]
width=max(150.,1.3*max(np.ptp(s['x'][ids]*1e6,axis=0).max() for ids in ids_all))
cols=2;rows=(len(labels)+1)//2;fig=plt.figure(figsize=(13,6*rows+1));fig.subplots_adjust(wspace=.18,hspace=.25,top=.87,bottom=.06)
for j,ids in enumerate(ids_all):
    xyz=s['x'][ids]*1e6;center=(xyz.min(0)+xyz.max(0))/2;internal=np.isin(b,ids).all(1)&(s['d']<1);cid=lookup[tuple((ids+1).tolist())]
    ax=ax3(fig,(rows,cols,j+1),f'{len(ids)} original nodes · {internal.sum()} live internal bonds\nComponent ID {cid}',(center-width/2,center+width/2))
    ax.scatter(*xyz.T,s=30,c='#37779b',depthshade=False);ax.add_collection3d(Line3DCollection(s['x'][b[internal]]*1e6,colors='#244e63',linewidths=.8,alpha=.5))
    records.append(dict(component_id=cid,original_atom_ids=ids+1,nodes=len(ids),internal_alive_bonds=int(internal.sum()),reference_shape=shape(ref['X'][ids]),actual_shape=shape(s['x'][ids])))
fig.suptitle(f'Completed uniform control · {float(s["time"])*1e3:.5f} ms · {len(labels)} true target-band free components\nSeparate local camera views; actual coordinates ×1; equal displayed width {width+2*max(.06*width,20.):.1f} µm\nBiphasic comparison pending — these are discrete connected point clouds, not reconstructed material shells',fontsize=15)
out=F/'previews';out.mkdir(exist_ok=True,parents=True)
for ext in ['png','pdf']:fig.savefig(out/('U_COMPLETE_TARGET_DETAILS.'+ext),dpi=250,bbox_inches='tight')
plt.close(fig)
dump(O/'U_COMPLETE_TARGET_DETAILS_EVIDENCE.json',dict(case='uniform',time_s=float(s['time']),scope='Completed uniform control only; no biphasic inference.',common_displayed_width_um=width+2*max(.06*width,20.),records=records))
print('UNIFORM_CONTROL_PREVIEW_DONE',flush=True)
