"""Read-only reconstruction/positive-energy diagnostics; no phase evolution."""
import numpy as np
from ..svk_spectral_split import green_lagrange, svk_split_energy


def bias_samples(discrete_F, analytic_F, material):
    """Ftilde callers must provide analytic F at bond quadrature midpoints."""
    psid=svk_split_energy(green_lagrange(discrete_F),material.lame,material.mu)[0]
    psia=svk_split_energy(green_lagrange(analytic_F),material.lame,material.mu)[0]
    return dict(F_error=np.linalg.norm(discrete_F-analytic_F,axis=(-2,-1)),
                Psi_plus_discrete=psid,Psi_plus_analytic=psia,Psi_plus_error=psid-psia)


def bias_summary(samples, weights, gradient_scale, psi_scale):
    e=samples['F_error'];p=samples['Psi_plus_error'];w=np.asarray(weights)/np.sum(weights)
    return dict(F_error_max=float(e.max()),F_error_max_scaled=float(e.max()/gradient_scale),
                F_error_RMS=float(np.sqrt(np.dot(w,e*e))),F_error_RMS_scaled=float(np.sqrt(np.dot(w,e*e))/gradient_scale),
                Psi_plus_error_max_Pa=float(abs(p).max()),Psi_plus_error_max_scaled=float(abs(p).max()/psi_scale),
                Psi_plus_error_RMS_Pa=float(np.sqrt(np.dot(w,p*p))),Psi_plus_error_RMS_scaled=float(np.sqrt(np.dot(w,p*p))/psi_scale),
                Psi_plus_error_mean_Pa=float(np.dot(w,p)))
