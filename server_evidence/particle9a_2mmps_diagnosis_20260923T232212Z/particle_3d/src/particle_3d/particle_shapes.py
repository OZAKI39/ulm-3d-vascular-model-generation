"""SI convex support geometry only; no flow, contact force, or membrane model."""
from dataclasses import dataclass, replace
import numpy as np
from .microbubble import vector3, positive_scalar
from .rbc_orientation import rotation_matrix

EPS = np.finfo(np.float64).eps
ROUNDOFF_FACTOR = 512


def unit(value):
    v = np.asarray(value, dtype=np.float64)
    if v.shape != (3,) or not np.isfinite(v).all():
        raise ValueError("finite shape (3,) direction required")
    length = np.linalg.norm(v)
    if length == 0 or not np.isfinite(length):
        raise ValueError("nonzero finite direction required")
    return v / length


def roundoff_length(*values):
    scale = max(float(np.max(np.abs(v))) for v in values)
    return ROUNDOFF_FACTOR * EPS * max(scale, np.finfo(float).tiny)


@dataclass(frozen=True)
class Sphere:
    center_m: np.ndarray
    radius_m: float
    mode: str = "SPHERE_MB"

    def __post_init__(self):
        object.__setattr__(self,"center_m",vector3(self.center_m,"center_m"))
        object.__setattr__(self,"radius_m",positive_scalar(self.radius_m,"radius_m"))

    @property
    def bounding_radius_m(self): return self.radius_m

    def support(self,direction): return self.center_m+self.radius_m*unit(direction)

    def moved(self,center_m): return replace(self,center_m=center_m)


@dataclass(frozen=True)
class Ellipsoid:
    center_m: np.ndarray
    axes_m: np.ndarray
    rotation: np.ndarray
    mode: str = "FREE_OBLATE"

    def __post_init__(self):
        object.__setattr__(self,"center_m",vector3(self.center_m,"center_m"))
        axes=vector3(self.axes_m,"axes_m")
        if np.any(axes<=0): raise ValueError("positive axes required")
        r=np.asarray(self.rotation,dtype=float)
        if r.shape!=(3,3) or not np.isfinite(r).all() or not np.allclose(r.T@r,np.eye(3),rtol=0,atol=256*EPS) or np.linalg.det(r)<0:
            raise ValueError("proper body-to-world rotation required")
        object.__setattr__(self,"axes_m",axes)
        object.__setattr__(self,"rotation",np.frombuffer(r.tobytes(),dtype=float).reshape(3,3))

    @classmethod
    def from_rbc(cls,center_m,geometry,q):
        return cls(center_m,[geometry.a_m,geometry.b_m,geometry.c_m],rotation_matrix(q))

    @property
    def bounding_radius_m(self): return float(np.max(self.axes_m))

    @property
    def quadratic(self): return (self.rotation*self.axes_m**2)@self.rotation.T

    def support(self,direction):
        d=self.rotation.T@unit(direction)
        return self.center_m+self.rotation@(self.axes_m**2*d/np.linalg.norm(self.axes_m*d))

    def moved(self,center_m): return replace(self,center_m=center_m)


@dataclass(frozen=True)
class Capsule:
    center_m: np.ndarray
    axis_world: np.ndarray
    radius_m: float
    cylindrical_length_m: float
    mode: str = "CAPILLARY_DEFORMED"

    def __post_init__(self):
        object.__setattr__(self,"center_m",vector3(self.center_m,"center_m"))
        object.__setattr__(self,"axis_world",vector3(unit(self.axis_world),"axis_world"))
        object.__setattr__(self,"radius_m",positive_scalar(self.radius_m,"radius_m"))
        length=float(self.cylindrical_length_m)
        if not np.isfinite(length) or length<0: raise ValueError("finite L>=0 required")
        object.__setattr__(self,"cylindrical_length_m",length)

    @property
    def bounding_radius_m(self): return self.radius_m+self.cylindrical_length_m/2

    @property
    def volume_m3(self): return np.pi*self.radius_m**2*self.cylindrical_length_m+4*np.pi*self.radius_m**3/3

    @property
    def area_m2(self): return 2*np.pi*self.radius_m*self.cylindrical_length_m+4*np.pi*self.radius_m**2

    def support(self,direction):
        d=unit(direction)
        return self.center_m+np.sign(d@self.axis_world)*self.cylindrical_length_m/2*self.axis_world+self.radius_m*d

    def moved(self,center_m): return replace(self,center_m=center_m)
