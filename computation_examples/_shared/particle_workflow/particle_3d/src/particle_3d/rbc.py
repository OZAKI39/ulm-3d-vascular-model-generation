"""Immutable rigid oblate RBC geometry/state; all core dynamics use SI.

q=(w,x,y,z), R(q): body→world; the body short axis is e3, so p=R(q)@e3.
"""
from dataclasses import dataclass
import numpy as np
from .microbubble import vector3, positive_scalar
from .rbc_distribution import DISTRIBUTION_ID
from .rbc_orientation import normalize_quaternion, short_axis


@dataclass(frozen=True, slots=True)
class GeometryProvenance:
    distribution_id: str
    contract_sha256: str
    seed: int
    rbc_id: int
    candidate_id: int
    D_um: float
    V_fL: float
    role: str


@dataclass(frozen=True, slots=True)
class RBCGeometry:
    a_m: np.float64
    b_m: np.float64
    c_m: np.float64
    volume_m3: np.float64
    provenance: GeometryProvenance

    def __post_init__(self):
        for name in ["a_m", "b_m", "c_m", "volume_m3"]:
            object.__setattr__(self, name, positive_scalar(getattr(self, name), name))
        tol = 256*np.finfo(float).eps
        if not np.isclose(self.a_m, self.b_m, rtol=tol, atol=0) or self.c_m >= self.a_m:
            raise ValueError("RBC geometry requires a=b>c>0")
        reconstructed = 4*np.pi*self.a_m*self.b_m*self.c_m/3
        if not np.isclose(reconstructed, self.volume_m3, rtol=tol, atol=0):
            raise ValueError("RBC axes and volume disagree")
        p = self.provenance
        if not isinstance(p, GeometryProvenance):
            raise ValueError("immutable geometry provenance required")
        if p.role != "SYNTHETIC_ONLY":
            if p.distribution_id != DISTRIBUTION_ID or len(p.contract_sha256) != 64 or min(p.seed,p.rbc_id,p.candidate_id) < 0:
                raise ValueError("biological RBC must carry sampler provenance")
            if not (np.isclose(self.a_m, p.D_um*.5e-6, rtol=tol, atol=0) and
                    np.isclose(self.volume_m3, p.V_fL*1e-18, rtol=tol, atol=0)):
                raise ValueError("geometry differs from sampled D/V provenance")

    @property
    def r(self):
        return np.float64(self.c_m/self.a_m)

    @property
    def jeffery_lambda(self):
        return np.float64((self.r*self.r-1)/(self.r*self.r+1))

    @property
    def distribution_id(self):
        return self.provenance.distribution_id

    @classmethod
    def from_population(cls, population, index):
        row, meta = population.samples[index], population.metadata
        p = GeometryProvenance(meta["contract_name"], meta["contract_sha256"], meta["seed"],
                               int(row["rbc_id"]), int(row["candidate_id"]), float(row["D_um"]), float(row["V_fL"]), meta["role"])
        return cls(row["a_m"], row["b_m"], row["c_m"], row["volume_m3"], p)

    @classmethod
    def synthetic_only(cls, a_m, c_m):
        return cls(a_m, a_m, c_m, 4*np.pi*a_m*a_m*c_m/3,
                   GeometryProvenance("SYNTHETIC_ONLY", "", 0, 0, 0, 0., 0., "SYNTHETIC_ONLY"))


@dataclass(frozen=True, slots=True)
class RBCState:
    particle_id: int
    position_m: np.ndarray
    quaternion_wxyz: np.ndarray
    velocity_m_s: np.ndarray
    angular_velocity_s_inv: np.ndarray
    geometry: RBCGeometry

    def __post_init__(self):
        if isinstance(self.particle_id, (bool,np.bool_)) or not isinstance(self.particle_id,(int,np.integer)) or self.particle_id < 0:
            raise ValueError("particle_id must be a nonnegative integer")
        if not isinstance(self.geometry, RBCGeometry):
            raise ValueError("RBCGeometry required")
        for name in ["position_m","velocity_m_s","angular_velocity_s_inv"]:
            object.__setattr__(self, name, vector3(getattr(self,name),name))
        q = normalize_quaternion(self.quaternion_wxyz)
        if q.shape != (4,):
            raise ValueError("one RBC requires a shape (4,) quaternion")
        object.__setattr__(self, "quaternion_wxyz", np.frombuffer(q.tobytes(), dtype=np.float64))

    @property
    def short_axis(self):
        return vector3(short_axis(self.quaternion_wxyz), "short_axis")
