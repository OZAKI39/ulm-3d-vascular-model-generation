"""Accepted-step fluid accounting by original atom membership."""
from common import *

class FragmentFluidRecorder:
    def __init__(self,ref,surface):
        self.ref=ref;self.surface=surface;self.node_work=np.zeros((len(ref['X']),2));self.birth={}

    def advance(self,p,s):
        drift=s['x']-p['x']
        for j,key in enumerate(['external_drive_force','external_drag_force']):
            self.node_work[:,j]+=.5*np.sum((p[key]+s[key])*drift,axis=1)

    def prepare(self,z,x,v,labels):
        n=len(labels);nc=int(labels.max())+1;lab=np.r_[labels,np.full(len(x)-n,-1)]
        owner=self.surface.tri[:,0];flab=lab[owner];valid=flab>=0;active=z['active']&valid
        def binc(values,mask=valid):return np.bincount(flab[mask],weights=values[mask],minlength=nc)
        self.area=binc(z['area']*z['weight']);self.faces=binc(z['active'].astype(float))
        self.unresolved=binc(z['unresolved'].astype(float));self.candidate=binc(z['candidate'].astype(float))
        methods=self.surface.method;self.l2=binc((z['active']&(methods==1)).astype(float));self.l1=self.faces-self.l2
        slip=z['outer_velocity']-(z['face_velocity']-z['normal']*np.sum(z['face_velocity']*z['normal'],axis=1)[:,None])
        self.slip_numerator=binc(np.linalg.norm(slip,axis=1)*z['area']*z['weight'])
        self.cover=np.zeros(len(x),bool);self.cover[np.unique(z['W'][:,active].nonzero()[0])]=True
        self.force=z['force'];self.x=x;self.v=v

    def component(self,cid,ids):
        # All integration work follows original atoms, including ancestry.
        # A distinct exact member set begins its own additional lifetime counter.
        com=self.x[ids].mean(0);force=self.force[ids];total=self.node_work[ids].sum(0)
        if cid not in self.birth:self.birth[cid]=total.copy()
        label=self.labels[ids[0]];area=float(self.area[label]);unknown=int(self.unresolved[label]);candidate=int(self.candidate[label])
        return dict(fluid_effective_area_m2=area,fluid_active_elements=int(self.faces[label]),fluid_L1_facets=int(self.l1[label]),fluid_L2_patches=int(self.l2[label]),
            fluid_unresolved_candidates=unknown,fluid_candidate_elements=candidate,fluid_unresolved_candidate_fraction=unknown/candidate if candidate else None,
            fluid_loaded_node_fraction=float(self.cover[ids].mean()),fluid_force_N=force.sum(0),fluid_moment_COM_Nm=np.cross(self.x[ids]-com,force).sum(0),
            fluid_power_W=float(np.sum(force*self.v[ids])),fluid_area_mean_tangential_slip_m_s=float(self.slip_numerator[label]/area) if area>0 else None,
            member_integrated_drive_work_J=float(total[0]),member_integrated_drag_work_J=float(total[1]),
            work_since_exact_member_set_birth_J=float((total-self.birth[cid]).sum()),
            fluid_loading_status='FLUID_LOADING_UNRESOLVED' if area==0 else 'PARTIAL_RESOLVED_PLUS_UNRESOLVED' if unknown else 'RESOLVED_APPROXIMATE_QUADRATURE',
            surface_method='L1_AND_L2' if self.l1[label] and self.l2[label] else 'L2_DUAL_FACE' if self.l2[label] else 'L1_SUPPORTED_TRIANGLE' if self.l1[label] else 'UNRESOLVED')

    def set_state(self,z,x,v,labels):
        self.labels=labels;self.prepare(z,x,v,labels)

    def save(self,path):
        np.savez_compressed(path,node_work_drive_drag_J=self.node_work,component_ids=np.array(list(self.birth),int),
            member_set_birth_work_J=np.array(list(self.birth.values())).reshape(-1,2))
