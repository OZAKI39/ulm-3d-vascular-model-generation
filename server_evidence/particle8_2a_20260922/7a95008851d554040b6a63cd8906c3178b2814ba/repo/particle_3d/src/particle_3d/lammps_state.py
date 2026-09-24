"""P6 storage schema. Physical geometry remains owned by particle_3d."""
from dataclasses import dataclass,asdict
import numpy as np
from .particle_shapes import Sphere,Ellipsoid,Capsule,EPS
from .rbc_orientation import rotation_matrix

SCHEMA_VERSION='PARTICLE6_STATE_V1'
PARTICLE_TYPES={'MB':1,'RBC':2}
SHAPE_MODES={'SPHERE_MB':1,'FREE_OBLATE':2,'CAPILLARY_DEFORMED':3}
PROPERTIES=[('i_particle_type','type_code',1),('i_shape_mode','mode_code',1),
 ('d2_quat','q',4),('d2_velocity','velocity',3),('d2_omega','omega',3),('d_radius','radius',1),
 ('d2_axes','axes',3),('d2_capsule_axis','capsule_axis',3),('d_capsule_radius','capsule_radius',1),
 ('d_capsule_length','capsule_length',1),('d_bound_radius','bound_radius',1),('d2_rotation','rotation',9)]
PROPERTY_COMMAND='fix particle_state all property/atom '+' '.join(name+(f' {width}' if '2_' in name else '') for name,_,width in PROPERTIES)+' ghost yes'


@dataclass(frozen=True)
class BridgeParticle:
    particle_id:int
    type_code:int
    mode_code:int
    position:np.ndarray
    q:np.ndarray
    velocity:np.ndarray
    omega:np.ndarray
    radius:float
    axes:np.ndarray
    capsule_axis:np.ndarray
    capsule_radius:float
    capsule_length:float
    bound_radius:float
    rotation:np.ndarray

    def __post_init__(self):
        if not isinstance(self.particle_id,(int,np.integer)) or not 0<self.particle_id<2**31:
            raise ValueError('Positive stable 32-bit particle ID required')
        if self.mode_code not in SHAPE_MODES.values() or self.type_code!=(1 if self.mode_code==1 else 2):
            raise ValueError('INVALID_OR_INFEASIBLE_PARTICLE_TYPE_SHAPE')
        for key,n in [('position',3),('q',4),('velocity',3),('omega',3),('axes',3),('capsule_axis',3),('rotation',9)]:
            v=np.asarray(getattr(self,key),dtype=float)
            if v.shape!=(n,) or not np.isfinite(v).all():raise ValueError('Invalid '+key)
            object.__setattr__(self,key,np.frombuffer(v.tobytes(),dtype=float))
        if abs(np.linalg.norm(self.q)-1)>512*EPS:raise ValueError('Stored quaternion must be unit within float64 roundoff')
        scalars=[self.radius,self.capsule_radius,self.capsule_length,self.bound_radius]
        if not np.isfinite(scalars).all() or min(scalars)<0 or self.bound_radius<=0:raise ValueError('Invalid shape scalars')
        shape=self.shape()
        if self.bound_radius!=shape.bounding_radius_m:raise ValueError('Stored bounding radius inconsistent with physical geometry')

    @classmethod
    def from_shape(cls,particle_id,shape,*,q=None,velocity=None,omega=None):
        if shape.mode not in SHAPE_MODES:raise ValueError('DEFORMATION_SURROGATE_INFEASIBLE_CANNOT_ENTER_LAMMPS')
        z=np.zeros(3);mode=SHAPE_MODES[shape.mode]
        return cls(particle_id,1 if mode==1 else 2,mode,shape.center_m,
            np.array([1.,0,0,0]) if q is None else q,z if velocity is None else velocity,z if omega is None else omega,
            shape.radius_m if mode==1 else 0.,shape.axes_m if mode==2 else z,
            shape.axis_world if mode==3 else z,shape.radius_m if mode==3 else 0.,
            shape.cylindrical_length_m if mode==3 else 0.,shape.bounding_radius_m,
            shape.rotation.ravel() if mode==2 else np.zeros(9))

    def shape(self):
        if self.mode_code==1:return Sphere(self.position,self.radius)
        if self.mode_code==2:return Ellipsoid(self.position,self.axes,self.rotation.reshape(3,3))
        if abs(np.linalg.norm(self.capsule_axis)-1)>512*EPS:raise ValueError('Invalid capsule axis')
        shape=Capsule(self.position,self.capsule_axis,self.capsule_radius,self.capsule_length)
        # Same exact-axis preservation used by P4 pure translation. No second normalization.
        object.__setattr__(shape,'axis_world',self.capsule_axis)
        return shape

    def to_dict(self):return asdict(self)


def schema_contract():
    return dict(schema_version=SCHEMA_VERSION,particle_types=PARTICLE_TYPES,shape_modes=SHAPE_MODES,
        properties=[dict(lammps_name=n,state_field=f,columns=w) for n,f,w in PROPERTIES],fix_id='particle_state',
        property_command=PROPERTY_COMMAND,quaternion_convention='P2_HAMILTON_WXYZ_BODY_TO_WORLD',
        capsule_quaternion_role='RETAINED_QUATERNION_NOT_CAPSULE_ORIENTATION',
        capsule_geometry_authority=['capsule_axis','capsule_radius','capsule_length'],
        extra_rotation_matrix_role='P4_EXACT_ROTATION_STORAGE_WITHOUT_QUATERNION_ROUNDTRIP_RECONSTRUCTION',
        physical_geometry_authority='PARTICLE_3D',mass_role='LAMMPS_BRIDGE_PLACEHOLDER_ONLY',ghost_properties=True)
