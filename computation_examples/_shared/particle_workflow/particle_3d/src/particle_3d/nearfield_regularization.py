"""User-approved sphere normal near-field regularization, separate from P5.

CONTINUUM_HANDOFF_CONTACT marks this continuum model's lower admissible gap;
it is NOT atomistic/molecular solid contact, glycocalyx or a shell offset.
Raw geometry is never clipped. Only the lubrication denominator is capped.
"""
from dataclasses import dataclass
import numpy as np
from .hydrodynamic_resistance import PhysicalNearField,evaluate_block,require_sphere

CONTRACT_NAME='NEAR_FIELD_REGULARIZATION_V1'
CHI_FULL=.01
CHI_OFF=.05
H_MOLECULAR_FLOOR=2e-9
H_SUSPENSION_RATIO=1e-3
STATES=['FAR_FIELD','NEARFIELD_TRANSITION','FULL_LUBRICATION','CONTINUUM_HANDOFF_CONTACT',
        'GEOMETRIC_HARD_CONTACT','NONSPHERICAL_NOT_SUPPORTED','BELOW_HANDOFF_TRIAL_REJECTED']


@dataclass(frozen=True)
class NearFieldRegularizationV1:
    h_molecular_floor_m:float=H_MOLECULAR_FLOOR
    role:str='USER_APPROVED_V1'

    def __post_init__(self):
        if self.role=='USER_APPROVED_V1' and self.h_molecular_floor_m==H_MOLECULAR_FLOOR:return
        if self.role=='SENSITIVITY_ONLY' and self.h_molecular_floor_m in [1.5e-9,2e-9,3e-9]:return
        raise ValueError('Only frozen V1 or explicitly approved sensitivity floors are allowed')

    @staticmethod
    def reference_length(a,b=None):
        require_sphere(a)
        if b is None:return a.radius_m
        require_sphere(b)
        return min(a.radius_m,b.radius_m)

    def lower_handoff_gap(self,a_ref):
        if not np.isfinite(a_ref) or a_ref<=0:raise ValueError('Positive finite reference length required')
        molecular=self.h_molecular_floor_m;suspension=H_SUSPENSION_RATIO*a_ref
        return dict(h_lower_m=max(molecular,suspension),h_molecular_component_m=molecular,
            h_suspension_component_m=suspension,which_component_dominates='MOLECULAR' if molecular>suspension else 'SUSPENSION' if suspension>molecular else 'EQUAL')

    @staticmethod
    def activation_weight(chi):
        if not np.isfinite(chi):raise ValueError('Finite dimensionless gap required')
        if chi<=CHI_FULL:return 1.
        if chi>=CHI_OFF:return 0.
        s=(chi-CHI_FULL)/(CHI_OFF-CHI_FULL)
        return (1-s)**2*(1+2*s) # exactly the requested cubic; stable near the outer endpoint

    @staticmethod
    def activation_derivative(chi):
        if not np.isfinite(chi):raise ValueError('Finite dimensionless gap required')
        if chi<=CHI_FULL or chi>=CHI_OFF:return 0.
        s=(chi-CHI_FULL)/(CHI_OFF-CHI_FULL)
        return 6*s*(s-1)/(CHI_OFF-CHI_FULL)

    @staticmethod
    def effective_gap(h_geom,h_lower):return max(h_geom,h_lower)

    @staticmethod
    def interaction_state(h_geom,a_ref,h_lower,roundoff=0.):
        if h_geom<=roundoff:return 'GEOMETRIC_HARD_CONTACT'
        if h_geom-h_lower < -roundoff:return 'BELOW_HANDOFF_TRIAL_REJECTED'
        if abs(h_geom-h_lower)<=roundoff:return 'CONTINUUM_HANDOFF_CONTACT'
        chi=h_geom/a_ref
        if chi>=CHI_OFF:return 'FAR_FIELD'
        if chi>CHI_FULL:return 'NEARFIELD_TRANSITION'
        return 'FULL_LUBRICATION'

    def evaluate(self,spec,particles,mu):
        if type(spec) is not PhysicalNearField:raise TypeError('Only exact-geometry PhysicalNearField specifications allowed')
        if not np.isfinite(mu) or mu<=0:raise ValueError('Positive finite viscosity required')
        a=particles[spec.particle_i_id];b=None if spec.particle_j_id is None else particles[spec.particle_j_id]
        if spec.particle_i_id==spec.particle_j_id:raise ValueError('Distinct pair IDs required')
        a_ref=self.reference_length(a,b);lower=self.lower_handoff_gap(a_ref);h=float(spec.gap_m);ro=float(spec.roundoff_m)
        if not np.isfinite([h,ro]).all() or ro<0:raise ValueError('Finite raw gap and nonnegative roundoff required')
        eff=self.effective_gap(h,lower['h_lower_m']);chi=h/a_ref;weight=self.activation_weight(chi)
        # Reuse the original P5 leading coefficient with its ORIGINAL R_eff.
        capped=PhysicalNearField(spec.particle_i_id,spec.particle_j_id,eff,spec.normal,0.,'VALIDATION_GAPS_ONLY')
        leading,n,old=evaluate_block(capped,particles,mu)
        coefficient=weight*leading;state=self.interaction_state(h,a_ref,lower['h_lower_m'],ro)
        record=dict(particle_i_id=spec.particle_i_id,particle_j_id=spec.particle_j_id,kind='WALL' if b is None else 'PAIR',
            model=CONTRACT_NAME,parameter_role=self.role,h_geom_m=h,gap_m=h,a_ref_m=a_ref,
            resistance_radius_m=old['relevant_radius_m'],chi=chi,w=weight,h_eff_m=eff,
            zeta_leading_kg_s=leading,coefficient_kg_s=coefficient,g_nf_m=h-lower['h_lower_m'],
            interaction_state=state,geometry_roundoff_m=ro,active=coefficient>0,
            continuum_state_admissible=h-lower['h_lower_m']>=-ro,normal=n.tolist(),**lower)
        return coefficient,n,record


def contract():
    return dict(contract_name=CONTRACT_NAME,status='USER_APPROVED_FOR_PARTICLE6_5',
        near_field_scope=['SPHERE_WALL','SPHERE_SPHERE'],chi_full=CHI_FULL,chi_off=CHI_OFF,
        smooth_transition='CUBIC_C1_SMOOTHSTEP',smoothstep_role='NUMERICAL_COUPLING_CHOICE',
        h_molecular_floor_m=H_MOLECULAR_FLOOR,h_molecular_floor_role='CONTINUUM_HANDOFF_SCALE',
        h_suspension_ratio=H_SUSPENSION_RATIO,h_lower_definition='max(h_molecular_floor_m, h_suspension_ratio*a_ref)',
        sphere_wall_reference_length='SPHERE_RADIUS',sphere_pair_reference_length='MIN_RADIUS',
        pair_reference_choice_role='CONSERVATIVE_REFERENCE_LENGTH_CHOICE',
        sphere_pair_resistance_radius='HARMONIC_REDUCED_RADIUS_REFF_UNCHANGED',
        h_eff_definition='max(h_geom,h_lower)',constraint_gap_definition='h_geom-h_lower',raw_geometry_gap_preserved=True,
        handoff_meaning='CONTINUUM_MODEL_BOUNDARY_NOT_MOLECULAR_SOLID_CONTACT',capped_lubrication_retained_at_handoff=True,
        contact_metric='UNCHANGED_P5_RESISTANCE_METRIC',physical_time_integrator='UNCHANGED_P3_BINARY_SUBDIVISION',
        position_projection=False,geometry_offset_m=0.,GLYCOCALYX_WALL_MODEL='NOT_FROZEN',
        SONOVUE_SHELL_THICKNESS_CONTEXT_NM=4.,shell_context_role='LITERATURE_CONTEXT_ONLY',sonovue_shell_radius_offset=False,
        nonspherical_lubrication='NOT_FROZEN',nonspherical_extension_forbidden=True,interaction_states=STATES,
        near_field_regularization_v1_frozen=True,production_neighbor_cutoff_frozen=False,production_neighbor_skin_frozen=False,
        production_particle_timestep_frozen=False,glycocalyx_model_frozen=False,nonspherical_lubrication_frozen=False,
        formulation_scope='SPHERE_NORMAL_NEAR_FIELD_REGULARIZATION_V1',full_production_hydrodynamics=False,
        sensitivity_floors_m=[1.5e-9,2e-9,3e-9],sensitivity_role='SENSITIVITY_ONLY',
        literature_roles=dict(A='DIRECT_SUSPENSION_SIMULATION_PRECEDENT',B='CROSS_DOMAIN_CONTINUUM_SCALE_EVIDENCE',
            C='CROSS_DOMAIN_NANOCONFINEMENT_EVIDENCE',D='LIPID_INTERFACE_SCALE_CONTEXT',E='SONOVUE_INTERFACE_SCALE_CONTEXT'),
        no_cfd_executed=True,particle7_started=False)
