"""Material connectivity uses surviving original bonds, never contact proximity."""
import numpy as np
from scipy.sparse import coo_matrix,csgraph
from scipy.sparse.linalg import eigsh
def topology(ref,d):
 n=int(ref['clot_count']);b=ref['bonds'];t=ref['bond_types'];live=d<1
 bb=b[live&(t==1)];matrix=coo_matrix((np.ones(2*len(bb)),(np.r_[bb[:,0],bb[:,1]],np.r_[bb[:,1],bb[:,0]])),shape=(n,n)).tocsr()
 nc,labels=csgraph.connected_components(matrix,directed=False);sizes=np.bincount(labels,minlength=nc);attached=np.zeros(nc,bool)
 it=b[live&(t==2)];material=it[it<n];attached[labels[material]]=True
 return labels,sizes,attached
def boundary(ref,nodes):
 mask=np.zeros(len(ref['X']),bool);mask[nodes]=True;b=ref['bonds']
 return np.flatnonzero(mask[b[:,0]]^mask[b[:,1]]),np.flatnonzero(mask[b[:,0]]&mask[b[:,1]])
def internal_components(ref,nodes,d):
 out,inside=boundary(ref,nodes);b=ref['bonds'];inside=inside[(ref['bond_types'][inside]==1)&(d[inside]<1)]
 idx=np.full(len(ref['X']),-1);idx[nodes]=np.arange(len(nodes));e=idx[b[inside]]
 return csgraph.connected_components(coo_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(len(nodes),len(nodes))),directed=False)[0]
def compactness(x,dx):
 span=np.ptp(x,axis=0);ev=np.linalg.eigvalsh((x-x.mean(0)).T@(x-x.mean(0))/len(x))
 return dict(min_m=x.min(0),max_m=x.max(0),span_m=span,voxel_bbox_fill=float(len(x)*dx**3/np.prod(span+dx)),radius_gyration_m=float(np.sqrt(ev.sum())),PCA_aspect_ratio=float(np.sqrt((ev[-1]+dx**2/12)/(ev[0]+dx**2/12))))
def proposals(ref,s):
 """Diagnostic graph partitions; no changes to the physical graph.

 Screening frozen before inspecting outputs: >=3 nodes, <=half material,
 reference bbox fill >=.2; >=80% internal original bonds survive;
 >=50% boundary integrity lost, greater than internal mean damage;
 and weighted boundary/total degree <=.15. These are search heuristics.
 """
 n=int(ref['clot_count']);b=ref['bonds'];bt=ref['bond_types'];d=s['d'];bb=b[bt==1];q=1-d[bt==1]
 proposed={};screen=[]
 for strength in [0.,.1,.25,.5,.75,.9,.99]:
  active=(d[bt==1]<1)&(q>=strength);e=bb[active]
  nc,lab=csgraph.connected_components(coo_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(n,n)),directed=False)
  sizes=np.bincount(lab)
  screen.append(dict(retained_integrity_threshold=strength,components=nc,sizes=sizes.tolist(),nonlargest_ge34=int(np.sum(sizes[np.arange(nc)!=np.argmax(sizes)]>=34)),nonlargest_ge114=int(np.sum(sizes[np.arange(nc)!=np.argmax(sizes)]>=114))))
  for j in np.flatnonzero((sizes>=3)&(sizes<=n//2)):
   nodes=tuple(np.flatnonzero(lab==j));proposed.setdefault(nodes,[]).append('integrity_threshold_'+str(strength))
 # Spectral sweep complements the threshold-component test: it can detect
 # bottlenecks held by undamaged bridges that no damage threshold removes.
 labels,sizes,attached=topology(ref,d);largest=int(np.argmax(sizes));nodes=np.flatnonzero(labels==largest)
 idx=np.full(n,-1);idx[nodes]=np.arange(len(nodes));keep=(idx[bb[:,0]]>=0)&(idx[bb[:,1]]>=0)&(q>0);edge=idx[bb[keep]];w=q[keep]
 W=coo_matrix((np.r_[w,w],(np.r_[edge[:,0],edge[:,1]],np.r_[edge[:,1],edge[:,0]])),shape=(len(nodes),len(nodes))).tocsr();L=csgraph.laplacian(W,normed=True)
 vals,vec=eigsh(L,k=4,which='SM',v0=np.linspace(1,2,len(nodes)),tol=1e-7)
 degree=np.asarray(W.sum(1)).ravel()
 for mode in range(1,4):
  for sign in [1.,-1.]:
   order=np.argsort(sign*vec[:,mode],kind='stable');rank=np.empty(len(nodes),int);rank[order]=np.arange(len(nodes));lo=np.minimum(rank[edge[:,0]],rank[edge[:,1]]);hi=np.maximum(rank[edge[:,0]],rank[edge[:,1]])
   delta=np.zeros(len(nodes)+1);np.add.at(delta,lo+1,w);np.add.at(delta,hi+1,-w);cut=np.cumsum(delta)[1:];vol=np.cumsum(degree[order]);conductance=cut/np.maximum(np.minimum(vol,vol[-1]-vol),1e-30)
   valid=np.arange(33,len(nodes)//2);ranked=valid[np.argsort(conductance[valid],kind='stable')]
   chosen=[]
   for k in ranked:
    if any(abs(int(k)-j)<20 for j in chosen):continue
    chosen.append(int(k));ids=tuple(sorted(nodes[order[:k+1]]));proposed.setdefault(ids,[]).append('spectral_sweep_mode_'+str(mode))
    if len(chosen)==3:break
 return proposed,screen
