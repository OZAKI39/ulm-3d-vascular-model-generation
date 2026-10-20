"""Frozen V3 quadrature adapted to the existing irregular Y clot and lumen.

The lumen is a fixed geometric constraint, not a newly solved fluid field.
Sidewall work is a separate external channel, never counted as fluid work.
"""
from common import *
from core_surface import CurrentSurface as CoreSurface
from pd_clot.lammps_fragmentation.native_y_wall_query import NativeYWallQuery
from pd_clot.lammps_fragmentation.native_y_wall_contact import CONTACT

class CurrentSurface(CoreSurface):
    def __init__(self,ref,mode,theta=None):
        super().__init__(ref,mode,theta)
        self.wall=NativeYWallQuery()
        self._lumen_cache=None;self._wall_cache=None
        boundary=load_npz(O/'BOUNDARY_CLASSIFICATION.npz')
        assert np.array_equal(boundary['tri'],self.base.tri)
        self.boundary_class=boundary['reference_class'];self.reference_wall_clearance=boundary['reference_wall_clearance_m']
        self.base.boundary_filter=self.boundary_filter

    def boundary_filter(self,x,d,centroid,original_candidate):
        exposed=original_candidate&(self.boundary_class==1)
        attachment=original_candidate&(self.boundary_class==2)
        free=self.topology(d)['free'];eligible=attachment&free[self.base.tri[:,0]]
        aperture=np.zeros(len(attachment));ids=np.flatnonzero(eligible)
        if len(ids):
            q=self.wall.query(centroid[ids]);aperture[ids]=q['clearance_m']-self.reference_wall_clearance[ids]
        opened=eligible&(aperture>1e-6)
        return exposed,opened,eligible,aperture,attachment&~eligible

    def geometry(self,x,d):
        base=super().geometry(x,d)
        if self._lumen_cache is not None and self._lumen_cache[0] is base:
            return self._lumen_cache[1]
        z=base.copy();ids=np.flatnonzero(z['active']);valid=np.zeros(len(z['active']),bool)
        if len(ids):
            points=z['centroid'][ids]
            # One-sided current lumen membership, without requiring an infinite
            # outward ray to remain in a bounded vessel.
            cells,_=self.wall.flow.locate(np.vstack([points,points+3e-6*z['normal'][ids]]))
            valid[ids]=(cells[:len(ids)]>=0)&(cells[len(ids):]>=0)
        rejected=z['active']&~valid
        for key in ['active','original','accessible','new_open']:
            z[key]=z[key]&valid
        z['weight']=np.where(valid,z['weight'],0.)
        z['occluded']=z['occluded']|rejected
        z['face_class']=z['face_class'].copy();z['face_class'][rejected]=4
        z['outside_lumen_rejected']=rejected
        self._lumen_cache=(base,z)
        return z

    def wall_geometry(self,x):
        xx=x[:self.n]
        if self._wall_cache is not None and np.array_equal(self._wall_cache[0],xx):
            return self._wall_cache[1]
        q=self.wall.query(xx)
        vector=xx-q['nearest_wall_point_m'];distance=np.linalg.norm(vector,axis=1)
        normal=vector/np.maximum(distance[:,None],1e-300)
        normal*=np.where(q['inside_lumen'],1.,-1.)[:,None]
        zero=distance<1e-15;normal[zero]=q['inward_normal'][zero]
        q['energy_normal']=normal
        q['overlap']=np.maximum(CONTACT['effective_radius_m']-q['clearance_m'],0.)
        # This short application has no transport handoff. If a point reaches
        # an open port, fail explicitly rather than treating it as a sidewall.
        for port in self.wall.ports.values():
            assert np.all((xx-port['center'])@port['outward']<0), 'OPEN_PORT_REACHED_REQUIRES_TRANSPORT_HANDOFF'
        self._wall_cache=(xx.copy(),q)
        return q

    def evaluate(self,x,v,d,t):
        z=super().evaluate(x,v,d,t);q=self.wall_geometry(x)
        vn=np.sum(v[:self.n]*q['energy_normal'],axis=1)
        damping=np.where(q['overlap']>0,CONTACT['damping_N_s_m']*np.maximum(-vn,0.),0.)
        force=np.zeros_like(x)
        force[:self.n]=(CONTACT['stiffness_N_m']*q['overlap']+damping)[:,None]*q['energy_normal']
        z.update(wall_force=force,wall_clearance=q['clearance_m'],
                 wall_U=float(.5*CONTACT['stiffness_N_m']*np.dot(q['overlap'],q['overlap'])),
                 wall_damping_power=float(np.dot(damping,vn)))
        return z

    def work(self,p,s,order=4):
        drift=s['x']-p['x']
        return np.array([.5*np.sum((p[key]+s[key])*drift) for key in
                         ['external_drive_force','external_drag_force','external_wall_force']])

