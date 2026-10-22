"""Single-particle wall-constrained replays, with one-way P0 free-flow queries."""
from dataclasses import dataclass,replace
import numpy as np
from .particle_shapes import Sphere,Ellipsoid,Capsule,EPS
from .rbc_orientation import advance_orientation,short_axis,angular_velocity
from .wall_gap import wall_gap,touching_contacts
from .wall_contact import constrain_contacts,ContactConstraintError
from .rbc_capillary_surrogate import select_rbc_shape,capsule_length,oblate_area,area_feasible_interval
from .physical_time_refinement import TrialNeedsSubdivision,swept_clearance_certificate,quasistatic_capsule_path_certificate


class SurrogateCaseStopped(RuntimeError):
    def __init__(self,record):self.record=record;super().__init__(str(record))


@dataclass(frozen=True)
class WallMotionState:
    time_s: float
    shape: object
    quaternion_wxyz: np.ndarray
    geometry: object | None
    gap: object
    free_velocity_m_s: np.ndarray
    corrected_velocity_m_s: np.ndarray
    free_omega_s_inv: np.ndarray
    tetra_id: int = -1
    boundary_event: str = 'ACTIVE'
    contact_normal_error_m_s: float = 0.
    deformation_search_evaluations: int = 0

    @property
    def shape_mode(self):return self.shape.mode

    def to_dict(self):
        p=short_axis(self.quaternion_wxyz)
        record=dict(time_s=float(self.time_s),center_m=self.shape.center_m.tolist(),q=self.quaternion_wxyz.tolist(),p=p.tolist(),
                    shape_mode=self.shape.mode,wall_gap_m=float(self.gap.gap_m),wall_triangle_id=int(self.gap.wall_triangle_id),
                    wall_feature=self.gap.wall_feature,wall_point_m=self.gap.wall_point_m.tolist(),
                    particle_point_m=self.gap.particle_point_m.tolist(),normal_inward=self.gap.normal_inward.tolist(),
                    contact_state=self.gap.state,roundoff_m=float(self.gap.roundoff_m),
                    free_velocity_m_s=self.free_velocity_m_s.tolist(),corrected_velocity_m_s=self.corrected_velocity_m_s.tolist(),
                    free_omega_s_inv=self.free_omega_s_inv.tolist(),tetra_id=int(self.tetra_id),boundary_event=self.boundary_event,
                    normal_constraint_error_m_s=float(self.contact_normal_error_m_s),
                    deformation_search_evaluations=self.deformation_search_evaluations)
        if isinstance(self.shape,Sphere):record['radius_m']=float(self.shape.radius_m)
        if self.geometry is not None:
            g=self.geometry;record.update(rbc_id=g.provenance.rbc_id,a_m=float(g.a_m),b_m=float(g.b_m),c_m=float(g.c_m),
                original_volume_m3=float(g.volume_m3),area_budget_m2=oblate_area(g.a_m,g.c_m),
                jeffery_orientation_interpreted=self.shape.mode=='FREE_OBLATE')
        if isinstance(self.shape,Capsule):
            record.update(R_cap_m=float(self.shape.radius_m),L_cap_m=float(self.shape.cylindrical_length_m),
                          capsule_axis_world=self.shape.axis_world.tolist(),capsule_volume_m3=float(self.shape.volume_m3),
                          capsule_area_m2=float(self.shape.area_m2),area_ratio=float(self.shape.area_m2/record['area_budget_m2']))
        return record


def flow_values(field,position,q,geometry,mode):
    sample=field.sample(position)
    if not sample.inside_lumen:return None
    omega=.5*sample.vorticity_s_inv if geometry is None else angular_velocity(short_axis(q),geometry.jeffery_lambda,sample.velocity_gradient_s_inv,sample.vorticity_s_inv)
    # In deformed mode this retained free diagnostic is NOT a capsule angular velocity.
    return np.asarray(sample.velocity_m_s),np.asarray(omega),int(sample.tetra_id)


def initial_wall_state(center,q,field,wall,*,radius=None,geometry=None,allow_deformation=True):
    values=flow_values(field,center,q,geometry,'FREE_OBLATE' if geometry else 'SPHERE_MB')
    if values is None:raise ValueError('P0 says initial center is outside lumen')
    v,omega,cell=values;count=0
    if geometry is None:
        shape=Sphere(center,radius);gap=wall_gap(shape,wall,inside_lumen=True)
    elif allow_deformation:
        decision=select_rbc_shape(center,geometry,q,v,wall)
        if decision.shape is None:raise SurrogateCaseStopped(dict(status=decision.status,time_s=0.,center_m=list(center),q=list(q),rbc_id=geometry.provenance.rbc_id,
            original_oblate_gap_m=decision.original_oblate_gap_m,reason=decision.reason,area_budget_m2=decision.area_budget_m2,volume_m3=float(geometry.volume_m3),accepted_state=False))
        shape,gap,count=decision.shape,decision.gap,decision.search_evaluations
    else:
        shape=Ellipsoid.from_rbc(center,geometry,q);gap=wall_gap(shape,wall,inside_lumen=True)
    if gap.state=='PENETRATING':raise ValueError('initial finite-size particle intersects WALL; no position projection')
    return WallMotionState(0.,shape,np.asarray(q,dtype=float),geometry,gap,v,v.copy(),omega,cell,deformation_search_evaluations=count)


def contact_trial(field,wall,*,allow_deformation=True,boundary_classifier=None,on_accept=None):
    def trial(old,dt):
        values=flow_values(field,old.shape.center_m,old.quaternion_wxyz,old.geometry,old.shape_mode)
        if values is None:raise TrialNeedsSubdivision('Old center has no P0 field')
        free,omega,cell=values
        contacts=touching_contacts(old.shape,wall) if old.gap.state=='TOUCHING' else []
        try:velocity,_,error=constrain_contacts(old.shape,free,omega,contacts)
        except ContactConstraintError as error:
            raise TrialNeedsSubdivision(str(error),old.gap.gap_m,old.gap.wall_triangle_id) from error
        center=old.shape.center_m+dt*velocity;event='ACTIVE'
        if boundary_classifier is not None:
            hit=boundary_classifier.first_event(old.shape.center_m,center)
            if hit is not None:
                if not hit.role.startswith('OUTLET_'):
                    raise TrialNeedsSubdivision('Center boundary event '+hit.role,old.gap.gap_m,old.gap.wall_triangle_id)
                dt*=hit.segment_fraction;center=hit.position_m;event=hit.role
        q=advance_orientation(old.quaternion_wxyz,omega,dt) if isinstance(old.shape,Ellipsoid) else old.quaternion_wxyz.copy()
        values=flow_values(field,center,q,old.geometry,old.shape_mode)
        if values is None:raise TrialNeedsSubdivision('Trial center outside P0 lumen',old.gap.gap_m,old.gap.wall_triangle_id)
        next_free,next_omega,next_cell=values;count=0
        if old.geometry is None:
            shape=old.shape.moved(center);gap=wall_gap(shape,wall,inside_lumen=True)
        elif allow_deformation:
            decision=select_rbc_shape(center,old.geometry,q,next_free,wall)
            if decision.shape is None:
                raise SurrogateCaseStopped(dict(status=decision.status,time_s=old.time_s,attempted_time_s=old.time_s+dt,
                    attempted_center_m=np.asarray(center).tolist(),last_accepted_state=old.to_dict(),reason=decision.reason,
                    original_oblate_gap_m=decision.original_oblate_gap_m,accepted_state=False))
            shape,gap,count=decision.shape,decision.gap,decision.search_evaluations
        else:
            shape=Ellipsoid.from_rbc(center,old.geometry,q);gap=wall_gap(shape,wall,inside_lumen=True)
        if gap.state=='PENETRATING':raise TrialNeedsSubdivision('Trial finite-size penetration',gap.gap_m,gap.wall_triangle_id)
        path_start,path_end=old.shape,shape;mode_proof=''
        if type(path_start) is not type(path_end):
            # A quasi-static geometry replacement is not a membrane trajectory.
            # Recompute a complete capsule passage path. Both replacement endpoint
            # geometries must be clear; no old fixed-oblate sweep bound is reused.
            if isinstance(shape,Capsule):
                path_start=Capsule(old.shape.center_m,free,shape.radius_m,shape.cylindrical_length_m)
                if wall_gap(path_start,wall).state=='PENETRATING':
                    raise TrialNeedsSubdivision('Mode change requires a safe recomputed capsule start')
            else:
                path_end=Capsule(center,next_free,old.shape.radius_m,old.shape.cylindrical_length_m)
                if wall_gap(path_end,wall).state=='PENETRATING':
                    raise TrialNeedsSubdivision('Mode change requires a safe capsule release endpoint')
            mode_proof='; QUASI_STATIC_REPLACEMENT_WITH_RECOMPUTED_CAPSULE_PATH'
        if isinstance(path_start,Capsule) and old.geometry is not None:
            clear,proof,triangle=quasistatic_capsule_path_certificate(path_start,path_end,wall,area_feasible_interval(old.geometry)[0])
        else:
            clear,proof,triangle=swept_clearance_certificate(path_start,path_end,wall,rotation_angle=np.linalg.norm(omega)*dt)
        if not clear:raise TrialNeedsSubdivision(proof,gap.gap_m,triangle)
        new=WallMotionState(old.time_s+dt,shape,q,old.geometry,gap,next_free,velocity,next_omega,next_cell,event,error,count)
        if on_accept is not None:on_accept(new)
        return new,proof+mode_proof
    return trial
