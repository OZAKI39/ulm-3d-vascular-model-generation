"""V3 conservative, supplementary material-dual-face surface quadrature.

The V2 supported-simplex path is retained verbatim through composition. Extra
patches belong only to genuinely free components, use actual nonuniform mapped reference cells,
and are excluded near existing L1 facets. They never change material topology.
This is a partial approximate surface, not a reconstructed continuum shell.
"""
from common import *
from scipy.sparse import csr_matrix,hstack,diags
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import ctypes,subprocess

from supported_surface import CurrentSurface as OldSurface
import mapped_helpers as mapped
OPEN=1e-6; FULL_OPEN=5e-6
CLASS_NAMES=['UNRESOLVED','ORIGINAL_EXPOSED','NEW_CRACK_ACCESSIBLE','OPEN_BUT_INACCESSIBLE','CLOSED_OR_OCCLUDED','DETACHED_EXPOSED']
MAX_REFERENCE_CONDITION=100.;MAX_CURRENT_CONDITION=100.
MAX_WEIGHT_L1=2.;MIN_WEIGHT=-.25

def smooth_open(a):
    q=np.clip((a-OPEN)/(FULL_OPEN-OPEN),0,1);return q*q*(3-2*q)

class CurrentSurface:
    def __init__(self,ref,mode,theta=None):
        self.ref=ref;self.base=OldSurface(ref,mode,theta);self.N=self.base.N;self.n=self.base.n
        self.field=self.base.field;self.label=LABEL;self.boundaries=self.base.boundaries;self.name=mode
        self.L=len(self.base.tri);self.M=6*self.n;self.loaded_ids=np.arange(self.n)
        self.owner=np.repeat(np.arange(self.n),6);self.axis=np.tile(np.repeat(np.arange(3),2),self.n)
        self.sign=np.tile([-1.,1.],3*self.n);self.offset=np.zeros((self.M,3));self.offset[np.arange(self.M),self.axis]=self.sign*DX/2
        self.reference_normal=np.vstack([self.base.reference_normal,self.offset/(DX/2)])
        self.reference_area=np.r_[self.base.reference_area,np.full(self.M,DX**2)]
        mapped.initialize(self)
        self.tri=np.vstack([self.base.tri,np.repeat(self.owner[:,None],3,axis=1)]) # L2 is NOT a triangle; indexing compatibility only.
        self.method=np.r_[np.zeros(self.L,np.uint8),np.ones(self.M,np.uint8)]
        self.material=np.r_[self.base.material,np.ones(self.M,bool)]
        self.original=np.r_[self.base.original,np.zeros(self.M,bool)]
        self.ever_loaded=np.zeros(self.L+self.M,bool);self.initial_area=None;self.cumulative_new_area=0.
        self._accepted=None;self._trial=None;self._geometry=None;self.topology_builds=0
        grid=ref['regular_X'];self.lattice=np.rint((grid-grid.min(0))/DX).astype(int)
        self.lookup={tuple(q):i for i,q in enumerate(self.lattice)}
        self.neighbor=np.array([self.lookup.get(tuple(self.lattice[i]+np.eye(3,dtype=int)[self.axis[k]]*int(self.sign[k])),-1) for k,i in enumerate(self.owner)])
        self.original[self.L:]=self.neighbor<0
        self.reference_tree=cKDTree(ref['X'])
        self.bridge=self._reference_crossings()
        b=ref['bonds'];self.bi=np.r_[b[:,0],b[:,1]];self.bj=np.r_[b[:,1],b[:,0]]
        self.br=(ref['X'][self.bj]-ref['X'][self.bi])/DX
        self.reference_weight=1/np.sum(self.br*self.br,axis=1)
        self._load_ray_kernel()

    def _load_ray_kernel(self):mapped.load_kernel(self)
    def _reference_crossings(self):return mapped.reference_crossings(self)

    def topology(self,d):
        failed=d>=1
        for value in [self._accepted,self._trial]:
            if value is not None and np.array_equal(value['failed'],failed):return value
        alive=~failed;b=self.ref['bonds'];edges=b[alive]
        graph=csr_matrix((np.ones(2*len(edges)),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(self.N,self.N))
        nc,lab=connected_components(graph,directed=False);attached=np.zeros(nc,bool);attached[lab[self.n:]]=True
        free=~attached[lab];neighbor=np.maximum(self.neighbor,0)
        boundary=free[self.owner]&((self.neighbor<0)|(lab[self.owner]!=lab[neighbor]))
        bridges=np.asarray(self.bridge@alive.astype(np.int32)).ravel()
        ii=self.bi[np.r_[alive,alive]];jj=self.bj[np.r_[alive,alive]];r=self.br[np.r_[alive,alive]];w=self.reference_weight[np.r_[alive,alive]]
        moment=np.zeros((self.N,3,3))
        for k in range(3):
            for l in range(3):moment[:,k,l]=np.bincount(ii,weights=w*r[:,k]*r[:,l],minlength=self.N)
        eig=np.linalg.eigvalsh(moment);good=(eig[:,0]>1e-10)&(eig[:,2]<=MAX_REFERENCE_CONDITION*eig[:,0])
        inverse=np.zeros_like(moment);inverse[good]=np.linalg.inv(moment[good])
        c=w[:,None]*np.einsum('nij,nj->ni',inverse[ii],r)
        gradients=[]
        for k in range(3):
            g=csr_matrix((c[:,k],(ii,jj)),shape=(self.N,self.N));g-=diags(np.asarray(g.sum(1)).ravel());gradients.append(g)
        # Face centroid weights reproduce affine position and its virtual velocity.
        # Rank-poor and strongly extrapolating patches remain unresolved.
        eligible=boundary&(bridges==0)&good[self.owner]
        faces=np.flatnonzero(eligible);rr=[];cc=[];vv=[];weight_ok=np.zeros(self.M,bool)
        for face in faces:
            i=self.owner[face];js=[i];vals=[1.]
            for a,g in enumerate(gradients):
                lo,hi=g.indptr[i:i+2];js.extend(g.indices[lo:hi]);vals.extend((self.offset[face,a]/DX*g.data[lo:hi]).tolist())
            js=np.asarray(js);vals=np.asarray(vals)
            merged={}
            for j,value in zip(js,vals):merged[int(j)]=merged.get(int(j),0.)+float(value)
            vals=np.array(list(merged.values()));js=np.array(list(merged))
            if vals.min(initial=0)<MIN_WEIGHT or abs(vals).sum()>MAX_WEIGHT_L1:continue
            assert abs(vals.sum()-1)<1e-12 and np.all(lab[js]==lab[i])
            weight_ok[face]=True;rr.extend(js.tolist());cc.extend([int(face)]*len(js));vv.extend(vals.tolist())
        W=csr_matrix((vv,(rr,cc)),shape=(self.N,self.M))
        result=dict(failed=failed.copy(),labels=lab,attached=attached,free=free,boundary=boundary,bridge_count=bridges,
                    gradients=gradients,reference_fit_valid=good,W=W,weight_ok=weight_ok)
        self._trial=result;self.topology_builds+=1
        return result

    def geometry(self,x,d):
        top=self.topology(d)
        if self._geometry is not None and self._geometry[1] is top and np.array_equal(self._geometry[0],x):return self._geometry[2]
        z1=self.base.geometry(x,d);owners=self.owner;M=self.M
        # Sparse reference gradients are rebuilt only for a changed failed-bond set.
        deformation=np.stack([g@x/DX for g in top['gradients']],axis=2)
        sv=np.linalg.svd(deformation,compute_uv=False)
        valid=top['reference_fit_valid']&(np.linalg.det(deformation)>1e-8)&(sv[:,2]>1e-6)&(sv[:,0]<=MAX_CURRENT_CONDITION*sv[:,2])
        inverse=np.zeros_like(deformation);inverse[valid]=np.linalg.inv(deformation[valid])
        center=np.asarray(top['W'].T@x)
        av=np.linalg.det(deformation[owners])[:,None]*np.einsum('nji,nj->ni',inverse[owners],self.reference_area_vectors)
        area=np.linalg.norm(av,axis=1);normal=av/np.maximum(area[:,None],1e-300)
        geometric=valid[owners]&np.isfinite(area)&(area>DX**2*1e-10)&top['weight_ok']
        # Reference footprint exclusion: retain all L1 faces, suppress entire dual
        # cells whose reference box may intersect an active L1 triangle. A sphere
        # enclosing both shapes is deliberately conservative (not area matching).
        covered=np.zeros(self.N,bool)
        ids=np.flatnonzero(z1['active'])
        if len(ids):
            p=self.ref['X'][self.base.tri[ids]];c=p.mean(1);rad=np.linalg.norm(p-c[:,None],axis=2).max(1)+self.cell_reference_radius.max()+1e-12
            for group in self.reference_tree.query_ball_point(c,rad):covered[group]=True
        overlap_guard=covered[owners]
        candidate=top['boundary'];eligible=candidate&geometric&(top['bridge_count']==0)&~overlap_guard
        # Exact mapped-cell normal separation across an original adjacent face.
        # An unresolvable opposite cell is a circumsphere obstacle, not fluid.
        nb=np.maximum(self.neighbor,0)
        corner_offset=np.einsum('nij,nkj->nki',deformation,self.cell_corner_offset)
        half_self=np.einsum('nkj,nj->nk',corner_offset[owners],normal).max(1)
        half_other=-np.einsum('nkj,nj->nk',corner_offset[nb],normal).min(1)
        half_other=np.where(valid[nb],half_other,self.cell_reference_radius[nb])
        aperture=np.sum((x[nb]-x[owners])*normal,axis=1)-half_self-half_other
        original=self.neighbor<0;aperture[original]=FULL_OPEN
        opened=eligible&(aperture>OPEN)
        radii=np.linalg.norm(corner_offset,axis=2).max(1);radii[~valid]=self.cell_reference_radius[~valid]
        ray_ids=np.flatnonzero(opened);closed=np.zeros(M,bool);blocked=np.zeros(M,bool)
        if len(ray_ids):
            inputs,cl,bl=mapped.mapped_ray(self,ray_ids,center,normal,owners,x,deformation,valid,radii)
            # Also respect the exact L1 solid support; the affine particle cells
            # are not assumed to enclose a nonaffinely distorted surviving tetra.
            a=self.base.owners[:,0];b=np.maximum(self.base.owners[:,1],0);solid=z1['supported_tetra']
            exterior=np.where(self.base.original,solid[a],solid[a]^solid[b])
            triangles=np.ascontiguousarray(x[self.base.tri[exterior]])
            self.raylib.ray_triangles(len(ray_ids),len(triangles),inputs[0].ctypes.data,inputs[1].ctypes.data,triangles.ctypes.data,3e-6,bl.ctypes.data)
            inside=self.base.occlusion_checker(x,solid)
            cl|=inside(inputs[0]+3e-6*inputs[1]).astype(np.uint8);bl|=cl
            closed[ray_ids]=cl;blocked[ray_ids]=bl
        accessible=opened&~blocked
        weight=np.where(accessible,np.where(original,1.,smooth_open(aperture)),0.)
        classes=np.zeros(M,np.uint8);classes[candidate&eligible&~opened]=4;classes[opened&blocked&~closed]=3;classes[closed]=4
        classes[accessible&original]=1;classes[accessible&~original]=5
        l1classes=np.zeros(self.L,np.uint8);l1classes[z1['original']]=1;l1classes[z1['accessible']]=2;l1classes[z1['new_open']&~z1['accessible']]=3;l1classes[z1['occluded']]=4
        z2=dict(centroid=center,normal=normal,area=area,weight=weight,active=weight>0,live=top['bridge_count']==0,geometric=geometric,
                quality=np.where(geometric,sv[owners,2]/np.maximum(sv[owners,0],1e-300),0.),area_ratio=area/self.reference_area[self.L:],
                original=accessible&original,new_open=opened&~original,accessible=accessible&~original,
                unresolved=candidate&(classes==0),candidate=candidate,aperture_m=aperture,fluid_region=np.where(accessible,0,-1),occluded=closed)
        out={k:np.concatenate([z1[k],z2[k]]) for k in z2}
        out.update(level1=z1,topology=top,face_class=np.r_[l1classes,classes],l2_overlap_guard=overlap_guard,
                   l2_bridge_count=top['bridge_count'],l2_cell_valid=valid,l2_deformation=deformation,
                   supported_tetra=z1['supported_tetra'],open_gap_tetra=z1['open_gap_tetra'],fluid_tetra=z1['fluid_tetra'],
                   gap_occluded=z1['gap_occluded'],unresolved_tetra=z1['unresolved_tetra'],volume_ratio=z1['volume_ratio'],
                   W=hstack([self.base.W,top['W']],format='csr'))
        self._geometry=(x.copy(),top,out);return out

    def evaluate(self,x,v,d,t):
        z=self.geometry(x,d).copy();W=z['W'];normal=z['normal'];vf=np.asarray(W.T@v)
        outer=self.field(z['centroid'],normal,t);vt=vf-normal*np.sum(vf*normal,axis=1)[:,None]
        td=BETA*outer;tr=-BETA*vt;area=z['area']*z['weight'];fd=area[:,None]*td;fr=area[:,None]*tr
        drive=np.asarray(W@fd);drag=np.asarray(W@fr);force=drive+drag
        z.update(face_velocity=vf,outer_velocity=outer,drive_traction=td,drag_traction=tr,traction=(td+tr)*z['weight'][:,None],face_force=fd+fr,
                 force=force,drive_force=drive,drag_force=drag,drive_power=float(np.sum(fd*vf)),drag_power=float(np.sum(fr*vf)),moment=np.cross(x,force).sum(0))
        return z

    def commit(self,z):
        self._accepted=z['topology'];new=z['active']&~self.ever_loaded
        self.cumulative_new_area+=float(z['area'][new].sum());self.ever_loaded|=z['active']
        if self.initial_area is None:self.initial_area=float((z['area']*z['weight']).sum())

    def rollback_trial(self):
        self._trial=None;self._geometry=None

    def save_checkpoint(self,path,state):
        assert self._accepted is not None and np.array_equal(self._accepted['failed'],state['d']>=1)
        np.savez_compressed(path,ever_loaded=self.ever_loaded,initial_area=np.array(self.initial_area),
            cumulative_new_area=np.array(self.cumulative_new_area),failed=self._accepted['failed'],
            step=np.array(state['step']),time=np.array(state['time']),atom_ids=state['atom_ids'],
            material_hash=np.array(sha(O/'provenance/reference.npz')),surface_version=np.array('V4_MAPPED_DUAL_FACE_1'))

    def restore_checkpoint(self,path,state):
        q=load_npz(path)
        assert str(q['material_hash'])==sha(O/'provenance/reference.npz') and str(q['surface_version'])=='V4_MAPPED_DUAL_FACE_1'
        assert int(q['step'])==state['step'] and abs(float(q['time'])-state['time'])<1e-14
        assert np.array_equal(q['atom_ids'],state['atom_ids']) and np.array_equal(q['failed'],state['d']>=1)
        assert q['ever_loaded'].shape==self.ever_loaded.shape
        self._accepted=self.topology(state['d']);self.ever_loaded=q['ever_loaded'].copy()
        self.initial_area=float(q['initial_area']);self.cumulative_new_area=float(q['cumulative_new_area'])

    def work(self,p,s,order=4):return self.base.work(p,s,order)
    def phase(self,t):return self.base.phase(t)

    def metrics(self,z,x):
        q={k:(v[:self.L] if isinstance(v,np.ndarray) and len(v.shape) and len(v)==self.L+self.M else v) for k,v in z.items()}
        out=self.base.metrics(q,x);active=z['active'];l2=active[self.L:]
        out.update(active_surface_faces=int(active.sum()),active_current_area_m2=float(z['area'][active].sum()),effective_loaded_area_m2=float(np.sum(z['area']*z['weight'])),
            original_exposed_faces=int(z['original'].sum()),new_open_faces=int(z['new_open'].sum()),fluid_accessible_new_faces=int(z['accessible'].sum()),
            closed_or_unresolved_faces=int(z['unresolved'].sum()),occluded_faces=int(z['occluded'].sum()),new_accessible_area_m2=float(z['area'][z['accessible']].sum()),
            actual_peak_slip_shear_Pa=float(np.linalg.norm(z['traction'],axis=1).max(initial=0)),actual_peak_drive_shear_Pa=float(np.linalg.norm(z['drive_traction'][active],axis=1).max(initial=0)),
            max_surface_speed_m_s=float(np.linalg.norm(z['face_velocity'][active],axis=1).max(initial=0)),max_outer_speed_m_s=float(np.linalg.norm(z['outer_velocity'][active],axis=1).max(initial=0)),
            max_node_force_N=float(np.linalg.norm(z['force'],axis=1).max(initial=0)),L2_active_patches=int(l2.sum()),L2_effective_area_m2=float((z['area']*z['weight'])[self.L:].sum()),
            L2_unresolved_patches=int(z['unresolved'][self.L:].sum()),L2_candidate_patches=int(z['candidate'][self.L:].sum()),surface_topology_builds=self.topology_builds)
        return out
