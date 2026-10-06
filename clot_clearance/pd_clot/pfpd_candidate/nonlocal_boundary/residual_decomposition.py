"""Force, moment and work accounting, retaining quadrature and assembly terms."""
import numpy as np
from .diagnostics import virtual_variations


def balances(X, arrays, continuum, force_ref, length, *, variations=None, moment_source=None):
    fd, bc, body, model = [arrays[k] for k in ('f_FD', 'f_BC_localized', 'f_body_local', 'f_model')]
    def accounting(vectors, surface, exact_body, scale):
        fi, fc, fb, fm = (v.sum(axis=0) for v in vectors)
        quadrature = fb-exact_body; lhs = fc-surface
        rhs = fm-quadrature-fi
        return dict(surface_local=surface, body_local_exact=exact_body, body_local_nodal=fb,
                    internal=fi, BC_localized=fc, model=fm, body_quadrature_error=quadrature,
                    BC_minus_surface=lhs, model_minus_quadrature_minus_internal=rhs,
                    closure=lhs-rhs, closure_scaled=float(np.linalg.norm(lhs-rhs)/scale),
                    formula='BC - surface = model - (nodal_body - exact_body) - internal',
                    simple_BC_minus_surface_minus_model=lhs-fm)
    force = accounting([fd, bc, body, model], np.asarray(continuum['surface_resultant_N']),
                       np.asarray(continuum['body_resultant_N']), force_ref)
    moment = accounting([np.cross(X, v) for v in [fd, bc, body, model]],
                        np.asarray(continuum['surface_moment_Nm']), np.asarray(continuum['body_moment_Nm']), force_ref*length)
    moment['convention'] = 'Moments about reference centroid X=0. Angular momentum conserves current-position moment; reference moment is not assumed zero.'
    if moment_source is not None:
        source = np.asarray(moment_source)
        moment['continuum_reference_moment_source'] = source
        rhs = moment['model_minus_quadrature_minus_internal']-source
        moment['model_minus_quadrature_minus_internal_minus_source'] = rhs
        moment['closure'] = moment['BC_minus_surface']-rhs
        moment['closure_scaled'] = float(np.linalg.norm(moment['closure'])/(force_ref*length))
        moment['formula'] += ' - continuum_reference_moment_source'
    work = []
    if variations is None:
        variations=[]
        for name,A,c in virtual_variations(length)+[('quadratic_z',None,None)]:
            if A is None:
                du=np.zeros_like(X);du[:,2]=X[:,2]**2/length
            else:du=X@A.T+c
            variations.append((name,du))
    for name,du in variations:
        wi, wc, wb, wm = (float(np.sum(v*du)) for v in [fd, bc, body, model])
        corrected = float(np.sum((fd+bc+body)*du))
        work.append(dict(variation=name, internal_J=wi, localized_J=wc, local_body_J=wb,
                         corrected_residual_J=corrected, bulk_model_J=wm,
                         closure_J=corrected-wm, closure_scaled=abs(corrected-wm)/(force_ref*max(length,float(np.linalg.norm(du,axis=1).max())))))
    return force, moment, work
