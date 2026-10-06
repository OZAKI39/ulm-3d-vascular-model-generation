"""Analytic homogeneous local solutions; explicit SVK material required."""
from dataclasses import dataclass
import numpy as np
from ..energy import constitutive
from .geometry import FACE_NAMES


@dataclass(frozen=True)
class AffineSurrogate:
    F: np.ndarray
    material: object
    translation: tuple = (0., 0., 0.)

    def position(self, X): return np.asarray(X)@self.F.T+np.asarray(self.translation)

    def displacement(self, X): return self.position(X)-X

    def response(self):
        W, P = constitutive(self.F, self.material)
        E = .5*(self.F.T@self.F-np.eye(3)); S = self.material.lame*np.trace(E)*np.eye(3)+2*self.material.mu*E
        return dict(F=self.F, E=E, S=S, P=P, W=float(W))

    def traction(self, normals): return np.asarray(normals)@self.response()['P'].T

    def exact_face_integrals(self, geometry, P=None):
        P = self.response()['P'] if P is None else np.asarray(P); rows = []
        for j, face in enumerate(FACE_NAMES):
            axis = j//2; normal = geometry.normals[j]
            area = geometry.volume/geometry.lengths[axis]
            midpoint = geometry.center+normal*geometry.lengths[axis]/2
            t = P@normal; tn = np.dot(t, normal)*normal; tt = t-tn
            rows.append(dict(face=face, normal=normal, centroid=midpoint, area_m2=area,
                             traction=t, traction_normal=tn, traction_tangential=tt,
                             force=area*t, normal_force=area*tn, tangential_force=area*tt,
                             moment=np.cross(midpoint-geometry.center, area*t)))
        return rows

    def surface_nodal_load(self, geometry, X, faces=FACE_NAMES):
        load = np.zeros_like(X)
        for q in geometry.face_quadrature(X):
            if q['face'] in faces:
                load[q['ids']] += q['area'][:, None]*self.traction(q['normal'])
        return load


def pressure_shear_traction(pressure, normals, tangential):
    """Future input representation only; no CFD provider or wall law."""
    n = np.asarray(normals); tau = np.asarray(tangential)
    if not np.allclose(np.sum(tau*n, axis=-1), 0., atol=1e-12):
        raise ValueError('Shear input must be tangential')
    return -np.asarray(pressure)[..., None]*n+tau
