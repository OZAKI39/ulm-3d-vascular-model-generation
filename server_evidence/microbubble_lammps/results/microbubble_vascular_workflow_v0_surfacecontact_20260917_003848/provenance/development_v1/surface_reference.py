"""Independent NumPy topology/contact oracle; never calls engine contact code."""
from pathlib import Path
import numpy as np,itertools,collections
from scipy.spatial import cKDTree

def closest_triangles(x,t):
    ab=t[:,1]-t[:,0];ac=t[:,2]-t[:,0];ap=x-t[:,0]
    aa=np.einsum('ij,ij->i',ab,ab);bb=np.einsum('ij,ij->i',ab,ac);cc=np.einsum('ij,ij->i',ac,ac);pa=np.einsum('ij,ij->i',ap,ab);pc=np.einsum('ij,ij->i',ap,ac);det=aa*cc-bb*bb
    v=(cc*pa-bb*pc)/det;w=(aa*pc-bb*pa)/det;q=t[:,0]+v[:,None]*ab+w[:,None]*ac;inside=(v>=0)&(w>=0)&(v+w<=1)
    bary=np.column_stack([1-v-w,v,w]);best=np.where(inside,np.linalg.norm(x-q,axis=1),np.inf)
    for k in range(3):
        d=t[:,(k+1)%3]-t[:,k];f=np.clip(np.einsum('ij,ij->i',x-t[:,k],d)/np.einsum('ij,ij->i',d,d),0,1);e=t[:,k]+f[:,None]*d;dist=np.linalg.norm(x-e,axis=1);take=dist<best;best[take]=dist[take];q[take]=e[take];bary[take]=0;bary[take,k]=1-f[take];bary[take,(k+1)%3]=f[take]
    return best,q,bary

class SurfaceReference:
    def __init__(self,vertices,faces,regions=None):
        raw=vertices[faces];self.vertices,inv=np.unique(raw.reshape(-1,3),axis=0,return_inverse=True);face=inv.reshape(-1,3)
        order=sorted(range(len(face)),key=lambda i:tuple(sorted(face[i])))
        self.original_ids=np.array(order);self.faces=face[order];self.t=self.vertices[self.faces];self.regions=np.zeros(len(face),int) if regions is None else np.asarray(regions)[order]
        cross=np.cross(self.t[:,1]-self.t[:,0],self.t[:,2]-self.t[:,0]);self.area=np.linalg.norm(cross,axis=1)/2;self.normal=cross/(2*self.area[:,None]);assert np.all(self.area>0)
        self.centers=self.t.mean(axis=1);self.rmax=np.linalg.norm(self.t-self.centers[:,None],axis=2).max();self.tree=cKDTree(self.centers)
        self.edges=collections.defaultdict(list);self.incident=collections.defaultdict(list);self.adj=[[] for _ in self.faces];self.dihedral={};self.bad_faces=set()
        for i,f in enumerate(self.faces):
            for v in f:self.incident[int(v)].append(i)
            for k in range(3):self.edges[tuple(sorted((int(f[k]),int(f[(k+1)%3]))))].append(i)
        for edge,inc in self.edges.items():
            if len(inc)>2:self.bad_faces.update(inc)
            if len(inc)==2:
                i,j=inc;ang=float(np.degrees(np.arccos(np.clip(self.normal[i]@self.normal[j],-1,1))));self.adj[i].append(j);self.adj[j].append(i);self.dihedral[tuple(sorted(inc))]=ang
    @classmethod
    def production(cls,stage):
        g=dict(np.load(Path(stage)/'geometry/GEOMETRY_ARRAYS.npz'));mask=g['classes']<=2;return cls(g['points'],g['faces'][mask],g['classes'][mask])
    def one_ring(self,face,vertex,theta):
        allowed=set(self.incident[vertex]);found={face};stack=[face]
        while stack:
            i=stack.pop()
            for j in self.adj[i]:
                if j in allowed and j not in found and self.dihedral[tuple(sorted((i,j)))]<=theta:found.add(j);stack.append(j)
        return found
    def query(self,x,tie_factor=1.,theta=15.):
        x=np.asarray(x,float);_,nearest=self.tree.query(x);d0=closest_triangles(x,self.t[[nearest]])[0][0]
        tol=tie_factor*64*np.finfo(float).eps*max(np.linalg.norm(x),1e-6)
        ix=np.array(sorted(self.tree.query_ball_point(x,d0+self.rmax+tol)),int);d,p,bary=closest_triangles(x,self.t[ix]);minimum=float(d.min());pick=d<=minimum+tol;ix=ix[pick];d=d[pick];p=p[pick];bary=bary[pick]
        candidates=[]
        for face,dist,pt,b in zip(ix,d,p,bary):
            face=int(face);zero=np.flatnonzero(b<=1e-9);support={face};feature='FACE';vertex=-1;edge=[]
            if len(zero)>=2:
                feature='VERTEX';vertex=int(self.faces[face,np.argmax(b)]);support=self.one_ring(face,vertex,theta)
            elif len(zero)==1:
                feature='EDGE';edge=sorted(int(self.faces[face,k]) for k in range(3) if k!=zero[0]);support={j for j in self.edges[tuple(edge)] if j==face or tuple(sorted((face,j))) in self.dihedral and self.dihedral[tuple(sorted((face,j)))]<=theta}
            candidates.append(dict(id=face,original_id=int(self.original_ids[face]),distance=float(dist),point=pt.tolist(),barycentric=b.tolist(),feature=feature,vertex=vertex,edge=edge,support=sorted(support),region=int(self.regions[face])))
        if any(c['id'] in self.bad_faces for c in candidates):return dict(status='FAIL_NONMANIFOLD',candidates=candidates,clusters=[],minimum_distance=minimum,tie_tolerance=tol)
        parent=list(range(len(ix)))
        def root(i):
            while parent[i]!=i:i=parent[i]
            return i
        for i,j in itertools.combinations(range(len(ix)),2):
            a,b=candidates[i],candidates[j];same=b['id'] in a['support'] or a['id'] in b['support']
            key=tuple(sorted((a['id'],b['id'])))
            if same or key in self.dihedral and self.dihedral[key]<=theta:parent[root(j)]=root(i)
        groups=collections.defaultdict(list)
        for i in range(len(ix)):groups[root(i)].append(i)
        clusters=[]
        for members in groups.values():
            supports=set();features=[];vertex_ids=set()
            for i in members:
                c=candidates[i];supports.update(c['support']);features.append(c['feature']);
                if c['vertex']>=0:vertex_ids.add(c['vertex'])
            vertex=next(iter(vertex_ids)) if len(vertex_ids)==1 and all(f=='VERTEX' for f in features) else -1
            faces=sorted(supports);weights=[]
            for f in faces:
                if vertex>=0:
                    k=list(self.faces[f]).index(vertex);a=self.t[f,(k+1)%3]-self.t[f,k];b=self.t[f,(k+2)%3]-self.t[f,k];weights.append(np.arccos(np.clip(a@b/(np.linalg.norm(a)*np.linalg.norm(b)),-1,1)))
                else:weights.append(self.area[f])
            normal=np.sum(np.array(weights)[:,None]*self.normal[faces],axis=0);length=np.linalg.norm(normal)
            if length<1e-14*sum(weights):raise ValueError('STOP_WALL_NORMAL_ORIENTATION_INVALID')
            normal/=length;closest=np.average(p[members],axis=0,weights=self.area[ix[members]])
            if normal@(x-closest)<0:normal=-normal
            if normal@(x-closest)<=0:raise ValueError('STOP_WALL_NORMAL_ORIENTATION_INVALID')
            mode='SMOOTH_VERTEX' if vertex>=0 else 'SMOOTH_EDGE' if 'EDGE' in features else 'SINGLE_FACE' if len(members)==1 else 'SINGLE_SURFACE_PATCH'
            clusters.append(dict(id=min(faces),candidate_ids=[int(ix[i]) for i in members],support=faces,normal=normal.tolist(),point=closest.tolist(),feature=mode))
        clusters.sort(key=lambda c:c['id'])
        return dict(status='TRUE_MULTI_SURFACE' if len(clusters)>1 else clusters[0]['feature'],candidates=candidates,clusters=clusters,minimum_distance=minimum,tie_tolerance=tol)

def project_velocity(raw,normals):
    raw=np.asarray(raw,float);N=np.asarray(normals,float).reshape(-1,3);scale=max(np.linalg.norm(raw),1e-30);tol=2e-12*scale
    best=None
    for rank in range(min(3,len(N))+1):
        for subset in itertools.combinations(range(len(N)),rank):
            if rank:
                A=N[list(subset)];G=A@A.T
                if np.linalg.matrix_rank(G,tol=1e-12)<rank:continue
                lam=np.linalg.solve(G,-A@raw)
                if min(lam)<-tol:continue
                used=raw+A.T@lam
            else:used=raw.copy();lam=np.array([])
            if len(N) and min(N@used)<-tol:continue
            objective=float(np.sum((used-raw)**2)/2)
            if best is None or objective<best[0]:best=(objective,used,subset,lam)
    if best is None:raise ValueError('STOP_MULTINORMAL_PROJECTION_FAILED')
    return best
