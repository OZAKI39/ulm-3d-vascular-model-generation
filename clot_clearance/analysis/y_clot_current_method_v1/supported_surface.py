"""Conservative supported-simplex surface / open-gap connectivity proxy.

No tetra is a material element: native mass, constitutive law and all original
bonds remain untouched. Lost support alone never turns a tetra into fluid.
Only locally disconnected, dilated, positive-volume cells can transmit fluid,
through dilated, disconnected triangular portals connected to the exterior.
New loaded triangles must bound a fully live positive-volume material tetra.
This is a lower-resolution support-complex proxy, NOT a continuum crack mesh.
"""
from common import *
from shear_field_library import ShearField
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from itertools import permutations, combinations

OPEN=1e-6
FULL_OPEN=5e-6
PAIR=np.array(list(combinations(range(4),2)))

def ramp(a):
    q=np.clip((a-OPEN)/(FULL_OPEN-OPEN),0.,1.)
    return q*q*(3-2*q)

class CurrentSurface:
    def __init__(self,ref,mode,theta=None):
        self.ref=ref;self.N=len(ref['X']);self.n=int(ref['clot_count']);self.name=mode
        X=ref['X']; axes=[np.unique(X[:,j]) for j in range(3)]
        lattice=np.rint((X-X.min(0))/DX).astype(int)
        lookup={tuple(q):i for i,q in enumerate(lattice)};tets=[]
        for key in lookup:
            origin=np.array(key)
            if any(key[j]>=len(axes[j])-1 for j in range(3)):continue
            for perm in permutations(range(3)):
                q=origin.copy();tet=[lookup[tuple(q)]]
                for j in perm:
                    q=q.copy();q[j]+=1;tet.append(lookup.get(tuple(q),-1))
                if min(tet)>=0:tets.append(tet)
        bonds={tuple(sorted(q)):i for i,q in enumerate(ref['bonds'])}
        # The historical Y graph omits some geometric neighbor pairs, especially
        # near the attachment. A support simplex requires all six ORIGINAL
        # edges; never manufacture an edge to complete a tetrahedron.
        self.tet=np.array([t for t in tets if all(tuple(sorted((t[i],t[j]))) in bonds for i,j in PAIR)],int);T=self.tet
        self.te=np.array([[bonds[tuple(sorted(t[p]))] for p in PAIR] for t in T])
        faces={}
        for ti,t in enumerate(T):
            for opposite in range(4):
                face=tuple(sorted(np.delete(t,opposite)))
                faces.setdefault(face,[]).append((ti,int(t[opposite])))
        self.tri=np.array(list(faces),int)
        self.owners=np.array([[q[0][0],q[1][0] if len(q)>1 else -1] for q in faces.values()])
        self.opposite=np.array([[q[0][1],q[1][1] if len(q)>1 else q[0][1]] for q in faces.values()])
        self.original=self.owners[:,1]<0
        self.material=np.all(self.tri<self.n,axis=1)
        self.edge_ids=np.array([[bonds[tuple(sorted((t[i],t[j])))] for i,j in [(0,1),(1,2),(2,0)]] for t in self.tri])
        self.surface_bonds=np.unique(np.r_[self.te.ravel(),self.edge_ids.ravel()])
        cross=np.cross(X[self.tri[:,1]]-X[self.tri[:,0]],X[self.tri[:,2]]-X[self.tri[:,0]])
        orient=np.sum(cross*(X[self.opposite[:,0]]-X[self.tri[:,0]]),axis=1)>0
        self.tri[orient]=self.tri[orient][:,[0,2,1]]
        cross=np.cross(X[self.tri[:,1]]-X[self.tri[:,0]],X[self.tri[:,2]]-X[self.tri[:,0]])
        self.reference_area=np.linalg.norm(cross,axis=1)/2
        self.reference_normal=cross/(2*self.reference_area[:,None])
        self.ref_volume=self.volume(X)
        self.ref_height=np.stack([np.abs(np.sum((X[self.opposite[:,j]]-X[self.tri[:,0]])*self.reference_normal,axis=1)) for j in range(2)],axis=1)
        self.W=csr_matrix((np.full(self.tri.size,1/3),(self.tri.ravel(),np.repeat(np.arange(len(self.tri)),3))),shape=(self.N,len(self.tri)))
        self.loaded_ids=np.arange(self.n)
        # Lookup local graph partitions for the 64 possible six-edge states.
        self.disconnected=np.zeros(64,bool)
        for mask in range(64):
            a=np.eye(4,dtype=bool)
            for k,(i,j) in enumerate(PAIR):a[i,j]=a[j,i]=not(mask&(1<<k))
            for k in range(4):a|=a[:,k,None]&a[None,k,:]
            self.disconnected[mask]=not a.all()
        self.field=ShearField(mode,theta,float(X[:self.n,2].max()))
        self.label=LABEL;self.boundaries=[RISE,HOLD,OFF,COHORT,END]
        self._cache=None;self.ever_loaded=np.zeros(len(self.tri),bool)
        self.initial_area=None;self.cumulative_new_area=0.

    def volume(self,x):
        p=x[self.tet];return np.einsum('ij,ij->i',np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),p[:,3]-p[:,0])/6

    def occlusion_checker(self,x,solid):
        """Actual current solid-tetra point containment; no reference-only visibility."""
        p=x[self.tet[solid]]
        if not len(p):return lambda points:np.zeros(len(points),bool)
        centers=p.mean(1);radius=np.linalg.norm(p-centers[:,None],axis=2).max()
        tree=cKDTree(centers)
        inverse=np.linalg.inv(np.transpose(p[:,1:]-p[:,:1],(0,2,1)))
        def inside(points):
            nearby=tree.query_ball_point(points,radius+1e-12)
            counts=np.array([len(a) for a in nearby]);pi=np.repeat(np.arange(len(points)),counts)
            if not len(pi):return np.zeros(len(points),bool)
            ti=np.concatenate(nearby).astype(int)
            bary=np.einsum('nij,nj->ni',inverse[ti],points[pi]-p[ti,0])
            hit=(bary>=-1e-9).all(1)&(bary.sum(1)<=1+1e-9)
            return np.bincount(pi,weights=hit,minlength=len(points))>0
        return inside

    def geometry(self,x,d):
        failed=d>=1
        if self._cache is not None and np.array_equal(x,self._cache[0]) and np.array_equal(failed,self._cache[1]):return self._cache[2].copy()
        vertices=x[self.tri];cross=np.cross(vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0])
        twice=np.linalg.norm(cross,axis=1);area=twice/2;n0=cross/np.maximum(twice[:,None],1e-300)
        ratio=area/self.reference_area
        lengths2=sum(np.sum((vertices[:,j]-vertices[:,(j+1)%3])**2,axis=1) for j in range(3))
        quality=4*np.sqrt(3)*area/np.maximum(lengths2,1e-300)
        geometric=(area>1e-12*DX**2)&np.isfinite(area)&(quality>1e-10)
        tfail=failed[self.te];mask=np.sum(tfail*(1<<np.arange(6)),axis=1)
        vr=self.volume(x)/self.ref_volume
        solid=(mask==0)&(vr>1e-8)
        gap_aperture=(vr-1)*DX
        gap=self.disconnected[mask]&(vr>0)&(gap_aperture>OPEN)
        inside=self.occlusion_checker(x,solid)
        gap_occluded=np.zeros(len(self.tet),bool)
        if gap.any():gap_occluded[gap]=inside(x[self.tet[gap]].mean(1))
        gap&=~gap_occluded
        ffailed=failed[self.edge_ids]
        # A triangle is disconnected iff at least two of its three edges failed.
        portal=(ffailed.sum(1)>=2)&geometric&((ratio-1)*DX>OPEN)
        if portal.any():portal[portal]&=~inside(vertices[portal].mean(1))
        a=self.owners[:,0];b=np.maximum(self.owners[:,1],0)
        internal=~self.original
        links=internal&portal&gap[a]&gap[b]
        graph=csr_matrix((np.ones(2*int(links.sum())),(np.r_[a[links],b[links]],np.r_[b[links],a[links]])),shape=(len(self.tet),len(self.tet)))
        _,regions=connected_components(graph,directed=False)
        seeds=self.original&portal&gap[a]
        fluid=gap&np.isin(regions,np.unique(regions[a[seeds]]))
        side=np.where(solid[a],0,1)
        owner=np.where(side==0,a,b);other=np.where(side==0,b,a)
        one_solid=np.where(self.original,solid[a],solid[a]^solid[b])
        normal=n0*np.where(side==0,1.,-1.)[:,None]
        face_live=~ffailed.any(1)
        normal_height=np.sum((x[self.opposite[np.arange(len(a)),1-side]]-vertices[:,0])*normal,axis=1)
        refheight=self.ref_height[np.arange(len(a)),1-side]
        aperture=normal_height-refheight
        original=self.original&solid[a]&face_live&geometric&self.material
        candidate=internal&one_solid&face_live&geometric&self.material
        opened=candidate&gap[other]&(aperture>OPEN)
        accessible=opened&fluid[other]
        visible=original|accessible
        occluded=np.zeros(len(self.tri),bool)
        if visible.any():
            points=vertices[visible].mean(1)
            # Two strictly outward current-space probes; neither may lie in solid.
            occluded[visible]=inside(points+1e-6*normal[visible])|inside(points+3e-6*normal[visible])
        original&=~occluded;accessible&=~occluded
        weight=np.where(original,1.,np.where(accessible,ramp(np.minimum(aperture,gap_aperture[other])),0.))
        active=weight>0
        unresolved=candidate&~accessible
        z=dict(centroid=vertices.mean(1),normal=normal,area=area,active=active,weight=weight,
            live=face_live,geometric=geometric,quality=quality,area_ratio=ratio,
            original=original,new_open=opened,accessible=accessible,unresolved=unresolved,
            candidate=candidate,aperture_m=aperture,fluid_region=np.where(accessible,regions[other],-1),
            occluded=occluded,gap_occluded=gap_occluded,
            supported_tetra=solid,open_gap_tetra=gap,fluid_tetra=fluid,
            unresolved_tetra=(mask!=0)&~fluid,volume_ratio=vr)
        self._cache=(x.copy(),failed.copy(),z)
        return z.copy()

    def evaluate(self,x,v,d,t):
        z=self.geometry(x,d);normal=z['normal'];vf=v[self.tri].mean(1)
        outer=self.field(z['centroid'],normal,t)
        vt=vf-normal*np.sum(vf*normal,axis=1)[:,None]
        td=BETA*outer;tr=-BETA*vt;area=z['area']*z['weight']
        fd=area[:,None]*td;fr=area[:,None]*tr
        drive=np.asarray(self.W@fd);drag=np.asarray(self.W@fr);force=drive+drag
        z.update(face_velocity=vf,outer_velocity=outer,drive_traction=td,drag_traction=tr,
            traction=(td+tr)*z['weight'][:,None],face_force=fd+fr,
            force=force,drive_force=drive,drag_force=drag,
            drive_power=float(np.sum(fd*vf)),drag_power=float(np.sum(fr*vf)),moment=np.cross(x,force).sum(0))
        return z

    def commit(self,z):
        new=z['active']&~self.ever_loaded
        self.cumulative_new_area+=float(z['area'][new].sum())
        self.ever_loaded|=z['active']
        if self.initial_area is None:self.initial_area=float(z['area'][z['active']].sum())

    def work(self,p,s,order=4):
        # Actual accepted endpoint forces dotted with actual native displacement.
        # Each channel is measured independently; energy residual is not hidden.
        dx=s['x']-p['x']
        return np.array([.5*np.sum((p[key]+s[key])*dx) for key in ['external_drive_force','external_drag_force']])

    def phase(self,t):return 'RISE' if t<RISE else 'HOLD' if t<HOLD else 'DRIVE_DOWN' if t<OFF else 'PASSIVE_DRAG_OBSERVATION'

    def metrics(self,z,x):
        active=z['active'];ratio=z['area_ratio'][active]
        original_den=float(self.reference_area[self.original&self.material].sum())
        return dict(active_surface_faces=int(active.sum()),inactive_broken_edge_faces=int((self.original&self.material&~z['live']).sum()),
            degenerate_live_faces=int((self.material&z['live']&~z['geometric']).sum()),
            active_reference_area_fraction=float(self.reference_area[z['original']].sum()/original_den),
            active_current_area_m2=float(z['area'][active].sum()),effective_loaded_area_m2=float(np.sum(z['area']*z['weight'])),
            original_exposed_faces=int(z['original'].sum()),new_open_faces=int(z['new_open'].sum()),fluid_accessible_new_faces=int(z['accessible'].sum()),
            closed_or_unresolved_faces=int(z['unresolved'].sum()),fluid_gap_cells=int(z['fluid_tetra'].sum()),
            occluded_faces=int(z['occluded'].sum()),occluded_gap_cells=int(z['gap_occluded'].sum()),
            unresolved_cells=int(z['unresolved_tetra'].sum()),inverted_cells=int((z['volume_ratio']<=0).sum()),
            new_accessible_area_m2=float(z['area'][z['accessible']].sum()),
            max_open_aperture_m=float(z['aperture_m'][z['new_open']].max(initial=0)),
            actual_peak_slip_shear_Pa=float(np.linalg.norm(z['traction'],axis=1).max(initial=0)),
            actual_peak_drive_shear_Pa=float(np.linalg.norm(z['drive_traction'][active],axis=1).max(initial=0)),
            max_active_area_ratio=float(ratio.max(initial=0)),min_active_area_ratio=float(ratio.min()) if len(ratio) else None,
            min_active_triangle_quality=float(z['quality'][active].min()) if active.any() else None,
            max_loaded_displacement_m=float(np.linalg.norm(x[:self.n]-self.ref['X'][:self.n],axis=1).max()),
            max_normal_rotation_deg=float(np.rad2deg(np.arccos(np.clip(np.sum(z['normal'][active]*self.reference_normal[active],axis=1),-1,1))).max(initial=0)),
            max_surface_speed_m_s=float(np.linalg.norm(z['face_velocity'][active],axis=1).max(initial=0)),
            max_outer_speed_m_s=float(np.linalg.norm(z['outer_velocity'][active],axis=1).max(initial=0)),
            max_node_force_N=float(np.linalg.norm(z['force'],axis=1).max(initial=0)))
