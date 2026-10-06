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


def stress_channel(P, channel='combined'):
    """Cartesian normal/tangential Piola projection, linear at frozen geometry."""
    P = np.asarray(P)
    if channel == 'combined': return P
    normal = P*np.eye(3)
    if channel == 'normal': return normal
    if channel == 'tangential': return P-normal
    raise ValueError(channel)


@dataclass(frozen=True)
class QuadraticShear:
    length: float = 18*.0000625
    gamma: float = .003
    beta_s: float = .001
    mu: float = 1000.

    @property
    def epsilon(self): return self.gamma  # Selective scaffold compatibility only.

    def displacement(self, X):
        X = np.asarray(X, float); u = np.zeros_like(X)
        u[..., 0] = self.gamma*X[..., 2]+self.beta_s*X[..., 2]**2/self.length
        return u

    def position(self, X): return np.asarray(X)+self.displacement(X)

    def response(self, Z):
        Z = np.asarray(Z, float); k = self.gamma+2*self.beta_s*Z/self.length; kp = 2*self.beta_s/self.length
        F = np.broadcast_to(np.eye(3), (*Z.shape, 3, 3)).copy(); F[..., 0, 2] = k
        E = np.zeros_like(F); E[..., 0, 2] = E[..., 2, 0] = k/2; E[..., 2, 2] = k*k/2
        P = np.zeros_like(F); P[..., 0, 0] = P[..., 2, 2] = self.mu*k*k
        P[..., 0, 2] = self.mu*(k+k**3); P[..., 2, 0] = self.mu*k
        b = np.zeros((*Z.shape, 3)); b[..., 0] = -self.mu*(1+3*k*k)*kp; b[..., 2] = -2*self.mu*k*kp
        return dict(F=F, C=np.eye(3)+2*E, E=E, S=2*self.mu*E, P=P,
                    W=self.mu/2*k*k+self.mu/4*k**4, b=b, J=np.ones_like(k), k=k, k_prime=kp)

    def body(self, Z, channel='combined'):
        b = self.response(Z)['b']
        if channel == 'normal': b[..., :2] = 0.
        elif channel == 'tangential': b[..., 2] = 0.
        elif channel != 'combined': raise ValueError(channel)
        return b

    def scales(self, h):
        r = self.response(np.array([-self.length/2, self.length/2]))
        p = np.linalg.norm(r['P'], axis=(-2, -1)).max(); b = np.linalg.norm(r['b'], axis=-1).max()
        return dict(P_ref_Pa=float(p), b_ref_N_m3=float(b), F_ref_N=float(max(p*h*h, b*h**3)),
                    gradient_scale=float(abs(self.gamma)+abs(self.beta_s)))

    def variations(self, X):
        X = np.asarray(X); result = []
        for name, A, c in virtual_variations(self.length):
            if name.startswith('translation_'): c = c/self.length  # Requested unit translations.
            result.append((name, X@A.T+c, np.broadcast_to(A, (len(X), 3, 3))))
        du = np.zeros_like(X); du[:, 0] = X[:, 2]**2/self.length
        grad = np.zeros((len(X), 3, 3)); grad[:, 0, 2] = 2*X[:, 2]/self.length
        result.append(('quadratic_x', du, grad))
        return result

    def continuum(self, channel='combined'):
        """Every face is integrated independently by polynomial-exact Gauss."""
        L = self.length; z, w = np.polynomial.legendre.leggauss(8); z *= L/2; w *= L/2
        ijk = np.indices((8,)*3).reshape(3, -1).T; X = z[ijk]; V = np.prod(w[ijk], axis=1)
        r = self.response(X[:, 2]); P = stress_channel(r['P'], channel); b = self.body(X[:, 2], channel)
        body = np.sum(V[:, None]*b, axis=0); mb = np.sum(np.cross(X, b)*V[:, None], axis=0)
        mb_current = np.sum(np.cross(self.position(X), b)*V[:, None], axis=0)
        moment_source = np.sum(np.cross(np.eye(3)[None], P.swapaxes(1, 2)).sum(axis=1)*V[:, None], axis=0)
        current_source = np.sum(np.cross(r['F'].swapaxes(1, 2), P.swapaxes(1, 2)).sum(axis=1)*V[:, None], axis=0)
        modes = self.variations(X); surface_work = np.zeros(len(modes)); faces = []
        ij = np.indices((8, 8)).reshape(2, -1).T; area = np.prod(w[ij], axis=1)
        for axis in range(3):
            for sign in [-1, 1]:
                Y = np.zeros((64, 3)); Y[:, axis] = sign*L/2; others = [a for a in range(3) if a != axis]
                Y[:, others] = z[ij]; N = sign*np.eye(3)[axis]
                t = stress_channel(self.response(Y[:, 2])['P'], channel)@N
                force = np.sum(t*area[:, None], axis=0)
                moment = np.sum(np.cross(Y, t)*area[:, None], axis=0)
                mc = np.sum(np.cross(self.position(Y), t)*area[:, None], axis=0)
                faces.append(dict(face='xyz'[axis]+('-' if sign<0 else '+'),normal=N,force_N=force,moment_Nm=moment,
                                  current_moment_Nm=mc,traction_min_Pa=t.min(axis=0),traction_max_Pa=t.max(axis=0)))
                for j, (_, du, _) in enumerate(self.variations(Y)): surface_work[j] += np.sum(t*du*area[:, None])
        surface = sum(f['force_N'] for f in faces); ms = sum(f['moment_Nm'] for f in faces)
        mc = sum(f['current_moment_Nm'] for f in faces)+mb_current
        sc = self.scales(L/18); rows = []
        for j, (name, du, grad) in enumerate(modes):
            bw = float(np.sum(b*du*V[:, None])); iw = float(np.sum(P*grad*V[:, None, None]))
            scale = sc['P_ref_Pa']*L*L*max(L,float(np.linalg.norm(du,axis=1).max()))
            rows.append(dict(variation=name,surface_J=surface_work[j],body_J=bw,stress_work_J=iw,
                             closure_J=surface_work[j]+bw-iw,closure_scaled=abs(surface_work[j]+bw-iw)/scale))
        exact_source = np.array([0.,self.mu*(self.gamma**3+self.gamma*self.beta_s**2)*L**3,0.])
        if channel=='normal': exact_source *= 0.
        return dict(channel=channel,faces=faces,surface_resultant_N=surface,body_resultant_N=body,
                    surface_moment_Nm=ms,body_moment_Nm=mb,reference_moment_source_Nm=moment_source,
                    reference_moment_source_analytic_Nm=exact_source,current_moment_Nm=mc,current_moment_source_Nm=current_source,
                    force_closure_scaled=float(np.linalg.norm(surface+body)/(sc['P_ref_Pa']*L**2)),
                    moment_closure_scaled=float(np.linalg.norm(ms+mb-exact_source)/(sc['P_ref_Pa']*L**3)),
                    moment_source_quadrature_error_scaled=float(np.linalg.norm(moment_source-exact_source)/(sc['P_ref_Pa']*L**3)),
                    current_moment_closure_scaled=float(np.linalg.norm(mc-current_source)/(sc['P_ref_Pa']*L**3)),
                    max_work_closure_scaled=max(v['closure_scaled'] for v in rows),work=rows,
                    energy_J=float(np.dot(V,r['W'])),integration='8-point tensor Gauss in volume and independently on each of six faces')
