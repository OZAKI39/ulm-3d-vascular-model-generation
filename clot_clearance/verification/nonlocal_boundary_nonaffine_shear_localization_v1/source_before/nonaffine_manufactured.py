"""Independent quadratic axial local reference; SI and reference coordinates.

No call to the general constitutive routine is used for these formulas.
"""
from dataclasses import dataclass
import numpy as np
from .diagnostics import virtual_variations


@dataclass(frozen=True)
class QuadraticAxial:
    length: float = .001125
    epsilon: float = .003
    beta: float = .001
    mu: float = 1000.

    def displacement(self, X):
        X = np.asarray(X); u = np.zeros_like(X)
        u[..., 2] = self.epsilon*X[..., 2]+self.beta*X[..., 2]**2/self.length
        return u

    def position(self, X): return np.asarray(X)+self.displacement(X)

    def response(self, Z):
        Z = np.asarray(Z); a = 1+self.epsilon+2*self.beta*Z/self.length
        if np.any(a <= 0): raise ValueError('Nonpositive reference Jacobian')
        F = np.broadcast_to(np.eye(3), (*Z.shape, 3, 3)).copy(); F[..., 2, 2] = a
        E = np.zeros_like(F); E[..., 2, 2] = .5*(a*a-1)
        S = 2*self.mu*E; P = np.zeros_like(F); P[..., 2, 2] = self.mu*(a**3-a)
        b = np.zeros((*Z.shape, 3)); b[..., 2] = -self.mu*(3*a*a-1)*(2*self.beta/self.length)
        return dict(F=F, E=E, S=S, P=P, W=self.mu/4*(a*a-1)**2, b=b, J=a)

    def scales(self, h):
        r = self.response(np.array([-self.length/2, self.length/2]))
        p = np.linalg.norm(r['P'], axis=(-2, -1)).max(); b = np.linalg.norm(r['b'], axis=-1).max()
        return dict(P_ref_Pa=float(p), b_ref_N_m3=float(b), F_ref_N=float(max(p*h*h, b*h**3)),
                    gradient_scale=float(abs(self.epsilon)+abs(self.beta)),
                    stress_gradient_scale_N_m3=float(2*self.mu*abs(2*self.beta/self.length)))

    def continuum(self):
        """8-point Gauss integrates every polynomial used here exactly."""
        L = self.length; z, w = np.polynomial.legendre.leggauss(8); z = z*L/2; w = w*L/2
        r = self.response(z); ends = self.response(np.array([-L/2, L/2]))
        surface = (ends['P'][1, :, 2]-ends['P'][0, :, 2])*L**2
        body = np.einsum('n,ni->i', w, r['b'])*L**2
        # Tensor quadrature also verifies all affine/shear/translation modes.
        ijk = np.indices((8,)*3).reshape(3, -1).T
        X = z[ijk]; V = np.prod(w[ijk], axis=1); rp = self.response(X[:, 2])
        rows = []
        for name, A, c in virtual_variations(L)+[('quadratic_z', None, None)]:
            if A is None:
                du = np.zeros_like(X); du[:, 2] = X[:, 2]**2/L
                grad = np.zeros((len(X), 3, 3)); grad[:, 2, 2] = 2*X[:, 2]/L
                end_du = np.array([0., 0., L/4]); sw = float(surface@end_du)
            else:
                du = X@A.T+c; grad = np.broadcast_to(A, (len(X), 3, 3))
                sw = sum(sign*L**2*ends['P'][k, :, 2]@(A@np.array([0., 0., sign*L/2])+c)
                         for k, sign in enumerate([-1, 1]))
            bw = float(np.sum(rp['b']*du*V[:, None])); iw = float(np.sum(rp['P']*grad*V[:, None, None]))
            rows.append(dict(variation=name, surface_J=float(sw), body_J=bw, stress_work_J=iw,
                             closure_J=float(sw+bw-iw)))
        scale = self.scales(L/18)
        return dict(surface_resultant_N=surface, body_resultant_N=body,
                    surface_moment_Nm=np.zeros(3), body_moment_Nm=np.zeros(3),
                    top_traction_Pa=ends['P'][1, :, 2], bottom_traction_Pa=-ends['P'][0, :, 2],
                    side_traction_Pa=np.zeros(3), energy_J=float(L**2*np.dot(w, r['W'])),
                    force_closure_N=float(np.linalg.norm(surface+body)),
                    force_closure_scaled=float(np.linalg.norm(surface+body)/(scale['P_ref_Pa']*L**2)),
                    work=rows, max_work_closure_scaled=max(abs(v['closure_J']) for v in rows)/(scale['P_ref_Pa']*L**3),
                    integration='8-point tensor Gauss; exact for the polynomial degrees present')
