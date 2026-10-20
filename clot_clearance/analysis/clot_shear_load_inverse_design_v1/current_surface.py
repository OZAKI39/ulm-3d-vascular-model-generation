"""Material-ID surface facets following current geometry and live edge support."""
from common import *
from shear_field_library import ShearField
from scipy.sparse import csr_matrix

class CurrentSurface:
    def __init__(self,ref,mode,theta=None):
        self.ref=ref;self.N=len(ref['X']);self.n=int(ref['clot_count']);self.name=mode
        data=load_npz(PREV/'SHEAR_TRACTION_FIELD.npz')
        self.tri=data['triangle_atom_ids'].astype(int)-1
        assert np.array_equal(data['reference_coordinates_m'],ref['X'])
        self.reference_area=data['area_m2'];self.reference_normal=data['n_clot']
        self.W=csr_matrix((np.full(self.tri.size,1/3),(self.tri.ravel(),np.repeat(np.arange(len(self.tri)),3))),shape=(self.N,len(self.tri)))
        edges=np.sort(np.stack([self.tri[:,[0,1]],self.tri[:,[1,2]],self.tri[:,[2,0]]],axis=1),axis=-1)
        lookup={tuple(sorted(pair)):i for i,pair in enumerate(ref['bonds'])}
        self.edge_ids=np.array([[lookup[tuple(e)] for e in face] for face in edges])
        self.surface_bonds=np.unique(self.edge_ids)
        assert not ref['initial_cut'][self.surface_bonds].any()
        assert (self.tri<self.n).all()
        self.field=ShearField(mode,theta,float(ref['X'][self.tri,2].mean()))
        self.label=LABEL;self.boundaries=[RISE,HOLD,OFF,COHORT,END]
        self.loaded_ids=np.unique(self.tri)

    def geometry(self,x,d):
        vertices=x[self.tri];cross=np.cross(vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0])
        twice=np.linalg.norm(cross,axis=1);area=twice/2
        normal=cross/np.maximum(twice[:,None],1e-300)
        lengths2=np.stack([np.sum((vertices[:,j]-vertices[:,(j+1)%3])**2,axis=1) for j in range(3)],axis=1)
        quality=4*np.sqrt(3)*area/np.maximum(lengths2.sum(1),1e-300)
        live=np.all(d[self.edge_ids]<1.,axis=1)
        geometric=(area>1e-12*DX**2)&np.isfinite(area)&(quality>1e-10)
        return dict(centroid=vertices.mean(1),normal=normal,area=area,active=live&geometric,
                    live=live,geometric=geometric,quality=quality,
                    area_ratio=area/self.reference_area)

    def evaluate(self,x,v,d,t):
        z=self.geometry(x,d);normal=z['normal'];vf=v[self.tri].mean(1)
        outer=self.field(z['centroid'],normal,t)
        vt=vf-normal*np.sum(vf*normal,axis=1)[:,None]
        td=BETA*outer;tr=-BETA*vt
        area=z['area']*z['active'];fd=area[:,None]*td;fr=area[:,None]*tr
        drive=np.asarray(self.W@fd);drag=np.asarray(self.W@fr);force=drive+drag
        z.update(face_velocity=vf,outer_velocity=outer,drive_traction=td,drag_traction=tr,
                 traction=(td+tr)*z['active'][:,None],face_force=fd+fr,
                 force=force,drive_force=drive,drag_force=drag,
                 drive_power=float(np.sum(fd*vf)),drag_power=float(np.sum(fr*vf)),
                 moment=np.cross(x,force).sum(0))
        return z

    def edge_damage_on_path(self,x,kprev):
        d=self.ref['initial_cut'].astype(float);ids=self.surface_bonds;b=self.ref['bonds'][ids]
        stretch=np.linalg.norm(x[b[:,1]]-x[b[:,0]],axis=1)/self.ref['r0'][ids]-1.
        k=np.maximum(kprev[ids],np.maximum(stretch,0.))
        d[ids]=(k>=SC[self.ref['bond_types'][ids]]).astype(float)
        return d

    def work(self,p,s,order=4):
        """Force dot actual native drift. Split exactly at facet-edge failure roots.

        The prescribed continuous interpolation is x=x0+a*(x1-x0), with
        v=(x1-x0)/dt. Irreversible edge history comes from the accepted old kappa.
        Native callbacks solve the linear slip force consistently with the
        final half-kick velocity without replacing the native integrator; quadrature
        disagreement remains a measured numerical diagnostic, never hidden work.
        """
        h=s['time']-p['time'];dx=s['x']-p['x'];v=dx/h
        ids=self.surface_bonds;new=ids[(p['d'][ids]<1)&(s['d'][ids]==1)]
        breaks=[0.,1.]
        if len(new):
            b=self.ref['bonds'][new];a=p['x'][b[:,1]]-p['x'][b[:,0]];c=dx[b[:,1]]-dx[b[:,0]]
            A2=np.sum(c*c,axis=1);B2=2*np.sum(a*c,axis=1)
            C2=np.sum(a*a,axis=1)-(self.ref['r0'][new]*(1+SC[self.ref['bond_types'][new]]))**2
            for aa,bb,cc in zip(A2,B2,C2):
                if aa>1e-50:
                    discriminant=np.sqrt(max(0.,bb*bb-4*aa*cc))
                    root=-2*cc/(bb+discriminant) if bb>=0 and bb+discriminant else (-bb+discriminant)/(2*aa)
                else:root=-cc/bb if bb else 1.
                if 1e-12<root<1-1e-12:breaks.append(float(root))
        breaks=np.unique(breaks);q,w=np.polynomial.legendre.leggauss(order);total=np.zeros(2)
        for lo,hi in zip(breaks[:-1],breaks[1:]):
            for qq,ww in zip(q,w):
                a=(lo+hi)/2+(hi-lo)*qq/2;x=p['x']+a*dx;t=p['time']+a*h
                d=self.edge_damage_on_path(x,p['kappa']);z=self.evaluate(x,v,d,t)
                total+=(hi-lo)*ww/2*np.array([np.sum(z['drive_force']*dx),np.sum(z['drag_force']*dx)])
        return total

    def phase(self,t):
        return 'RISE' if t<RISE else 'HOLD' if t<HOLD else 'DRIVE_DOWN' if t<OFF else 'PASSIVE_DRAG_OBSERVATION'

    def metrics(self,z,x):
        active=z['active'];ratio=z['area_ratio'][active];quality=z['quality'][active]
        return dict(active_surface_faces=int(active.sum()),inactive_broken_edge_faces=int((~z['live']).sum()),
                    degenerate_live_faces=int((z['live']&~z['geometric']).sum()),
                    active_reference_area_fraction=float(self.reference_area[active].sum()/self.reference_area.sum()),
                    active_current_area_m2=float(z['area'][active].sum()),
                    actual_peak_slip_shear_Pa=float(np.linalg.norm(z['traction'],axis=1).max(initial=0)),
                    actual_peak_drive_shear_Pa=float(np.linalg.norm(z['drive_traction'][active],axis=1).max(initial=0)),
                    max_active_area_ratio=float(ratio.max(initial=0)),min_active_area_ratio=float(ratio.min()) if len(ratio) else None,
                    min_active_triangle_quality=float(quality.min()) if len(quality) else None,
                    max_loaded_displacement_m=float(np.linalg.norm(x[self.loaded_ids]-self.ref['X'][self.loaded_ids],axis=1).max()),
                    max_normal_rotation_deg=float(np.degrees(np.arccos(np.clip(np.sum(z['normal'][active]*self.reference_normal[active],axis=1),-1,1))).max(initial=0)),
                    max_tangency_error_Pa=float(abs(np.sum(z['traction']*z['normal'],axis=1)).max(initial=0)))

def freeze_surface(ref):
    surface=CurrentSurface(ref,'uniform')
    np.savez_compressed(O/'ORIGINAL_SURFACE.npz',triangle_atom_ids=surface.tri+1,edge_original_bond_ids=surface.edge_ids+1,
                        reference_area_m2=surface.reference_area,reference_normal=surface.reference_normal)
    dump(O/'CURRENT_SURFACE_SPEC.json',dict(status='CURRENT_SURFACE_SHEAR_COUPLING',label=LABEL,
         original_surface='182 original outer top triangles, 112 original material vertices, inherited byte-identically from the prior sealed mapping. Other original faces are outside this model.',
         update='Every force callback uses current physical coordinates, current facet areas/normals, current material velocities and committed original-ID edge damage.',
         retention='A facet is retained only if all three original edge bonds remain live and its current area/quality is nondegenerate. This is conservative geometric support, not a fluid-accessibility proof.',
         detached='Attachment is not an input to the force law. Intact original exterior facets on detached components retain position-dependent drive and relative-velocity drag; outside the drive region only passive drag remains. No original dead traction survives by node ID.',
         unresolved='NEW_CRACK_FACE_WETTING_UNRESOLVED. Broken triangles are removed from the load discretization, never from the material graph. No new crack triangles, pressure, full fragment hydrodynamics or vessel transport are fabricated.',
         point_fragments='Fragments without any retained original surface triangle have unresolved fluid coupling; their force is not replaced by an invented area or automatic COM drag.',
         quadrature='Centroid traction and vertex-mean velocity, positive 1/3 nodal weights. Exact discrete resultant, current moment, and adjoint virtual power; not spatial mesh convergence.',
         numerical='Current callback solves the linear slip law consistently with the final native half-kick velocity, using a converged fixed point and unchanged native NVE. Continuous drift work uses actual dx/dt and exact edge-failure root splits; force-stage/end-velocity agreement and work errors are recorded.',
         diagnostics='Record first surface-support loss, displacement>0.1h, rotation>15deg, area ratio outside [.5,2], degenerate live faces, local Reynolds criteria. These flag separate limits, not a calibrated universal failure cutoff.',
         MPI_restart='Single rank only as inherited. Original IDs reorder callback arrays explicitly. External law is memoryless beyond native kappa/d and frozen parameters; native restart history/kernel unchanged. Only this callback path is checked.'))
