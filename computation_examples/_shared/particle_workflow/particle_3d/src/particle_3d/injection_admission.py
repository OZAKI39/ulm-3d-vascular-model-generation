"""Finite-body admission; open inlet has zero geometric offset."""
from collections import Counter
from copy import deepcopy
import numpy as np
from .particle_shapes import Sphere, Ellipsoid, roundoff_length
from .pair_geometry import pair_gap
from .wall_gap import wall_gap
from .rbc_capillary_surrogate import select_rbc_shape
from .injection_population import geometry_from_dict
from .nearfield_regularization import NearFieldRegularizationV1
from .lammps_state import BridgeParticle

MAX_CANDIDATE_DRAWS_PER_ADMISSION_ATTEMPT=512


def event_shape(event, position):
    if event['species']=='MB': return Sphere(position,event['radius_m'])
    return Ellipsoid.from_rbc(position,geometry_from_dict(event['geometry']),event['q'])


class FiniteSizeAdmission:
    def __init__(self,sampler,source,*,wall=None,field=None,velocity=(0.,0.,1e-3),guard=512):
        if not isinstance(guard,int) or not 1<=guard<=512: raise ValueError('Implementation safety guard must be in [1,512]')
        self.sampler=sampler; self.source=source; self.wall=wall; self.field=field; self.velocity=np.asarray(velocity,float)
        self.guard=guard; self.counts=Counter(); self.candidates=[]; self.policy=NearFieldRegularizationV1()
    def check(self,event,position,active):
        shape=event_shape(event,position); free=self.velocity; omega=np.zeros(3); detail={}
        if self.field is not None:
            sample=self.field.sample(position)
            if not sample.inside_lumen: return None,'CENTER_OUTSIDE_FROZEN_LUMEN',detail
            free=np.asarray(sample.velocity_m_s); omega=.5*np.asarray(sample.vorticity_s_inv)
            if event['species']=='RBC':
                from .rbc_orientation import angular_velocity,short_axis
                geometry=geometry_from_dict(event['geometry'])
                omega=angular_velocity(short_axis(event['q']),geometry.jeffery_lambda,sample.velocity_gradient_s_inv,sample.vorticity_s_inv)
                detail['omega_role']='P2_FREE_JEFFERY_DIAGNOSTIC; NOT_CAPSULE_PHYSICAL_ROTATION'
        if self.wall is not None:
            _,distance=self.wall.nearest_center_triangle(position)
            lower=self.policy.lower_handoff_gap(shape.radius_m)['h_lower_m'] if isinstance(shape,Sphere) else 0.
            clear_bound=distance-shape.bounding_radius_m
            if clear_bound>lower+self.wall.roundoff_m:
                # Certified separating bounding ball, not an effective RBC radius.
                # Avoid unnecessary narrow-phase solves on distant skinny triangles.
                detail['wall_gap_lower_bound_m']=float(clear_bound)
                detail['wall_certificate']='BOUNDING_BALL_SEPARATION'
            else:
                if isinstance(shape,Ellipsoid):
                    # P3 already tries the unchanged FREE_OBLATE first and has
                    # certified early conflict checks before its capsule search.
                    decision=select_rbc_shape(position,geometry_from_dict(event['geometry']),event['q'],free,self.wall)
                    detail['original_wall_gap_m']=float(decision.original_oblate_gap_m)
                    detail['deformation_status']=decision.status
                    if decision.shape is None: return None,'DEFORMATION_INFEASIBLE',detail
                    shape,gap=decision.shape,decision.gap
                else:
                    gap=wall_gap(shape,self.wall,inside_lumen=True)
                    detail['original_wall_gap_m']=float(gap.gap_m)
                if gap.state=='PENETRATING': return None,'WALL_REJECTED',detail
                detail['wall_gap_m']=float(gap.gap_m)
                if isinstance(shape,Sphere):
                    detail['wall_g_nf_m']=gap.gap_m-self.policy.lower_handoff_gap(shape.radius_m)['h_lower_m']
                    if detail['wall_g_nf_m'] < -gap.roundoff_m: return None,'WALL_NEARFIELD_REJECTED',detail
        # Bounding spheres only cull pairs with proven separation. Every potential
        # intersection is decided by the unchanged P4 exact convex pair-gap engine.
        for tag,other in active.items():
            other=other.shape() if isinstance(other,BridgeParticle) else other
            lower=self.policy.lower_handoff_gap(min(shape.radius_m,other.radius_m))['h_lower_m'] if isinstance(shape,Sphere) and isinstance(other,Sphere) else 0.
            ro=roundoff_length(shape.center_m,other.center_m,shape.bounding_radius_m,other.bounding_radius_m)
            if np.linalg.norm(shape.center_m-other.center_m)>shape.bounding_radius_m+other.bounding_radius_m+lower+ro: continue
            g=pair_gap(shape,other,event['particle_id'],tag)
            detail.update(pair_gap_m=float(g.gap_m),other_id=tag)
            if g.state=='PENETRATING': return None,'PAIR_REJECTED',detail
            if g.gap_m-lower < -g.roundoff_budget_m: return None,'PAIR_NEARFIELD_REJECTED',detail
        return BridgeParticle.from_shape(event['particle_id'],shape,q=event['q'],velocity=free,omega=omega),'ACCEPTED',detail
    def attempt(self,event,active):
        event['attempt_count']+=1
        name=event['species']+'_POSITION' if event['attempt_count']==1 else 'ADMISSION_RETRY'
        rng=self.source.rng[name]
        for draw in range(self.guard):
            position,triangle=self.sampler.sample(rng)
            p,status,detail=self.check(event,position[0],active)
            self.counts[status]+=1
            row=dict(particle_id=event['particle_id'],species=event['species'],attempt=event['attempt_count'],draw=draw,
                     position_m=position[0].tolist(),triangle_id=int(triangle[0]),status=status,**detail)
            self.candidates.append(row); event['last_candidate']=deepcopy(row)
            if p is not None: return p
        self.counts['GUARD_EXHAUSTED_PENDING']+=1
        return None


def plug_side_wall_margin(shape,width):
    """Exact support-plane lower gap for a rectangular common-translation fixture."""
    support=np.array([shape.support(n)[j]-shape.center_m[j] for j,n in enumerate(np.eye(3)[:2])])
    return float(np.min(width/2-np.abs(shape.center_m[:2])-support))


class PlugPathAdmission(FiniteSizeAdmission):
    """Synthetic-only admission that certifies the complete translated WALL path.

    Reject an uncertified path and try another inlet position with the same
    geometry. These are trajectory rejections, not measured WALL penetration.
    """
    def __init__(self,sampler,source,*,width,velocity,guard=512):
        super().__init__(sampler,source,velocity=velocity,guard=guard); self.width=width
    def check(self,event,position,active):
        shape=event_shape(event,position); margin=plug_side_wall_margin(shape,self.width)
        lower=self.policy.lower_handoff_gap(shape.radius_m)['h_lower_m'] if isinstance(shape,Sphere) else 0.
        if margin<lower+roundoff_length(self.width):
            return None,'FULL_PLUG_WALL_PATH_UNCERTIFIED',dict(wall_path_lower_bound_m=margin)
        p,status,detail=super().check(event,position,active)
        detail['wall_path_lower_bound_m']=margin
        return p,status,detail
