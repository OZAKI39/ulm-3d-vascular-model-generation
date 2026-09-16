"""Independent classical diagnostics, never inserted into package outputs.
Brenner series: Ascoli 1988 Caltech dissertation Eq.(29), PDF p23; reproduced
Brenner(1961). Published wall asymptotics: Sprinkle et al. 2020 Appendix A2,
Table I, arXiv:2005.06002 p37, explicitly distinguishes fitted linear terms.
"""
import mpmath as mp
import numpy as np

def normal_brenner(epsilon,dps=60):
    with mp.workdps(dps):
        e=mp.mpf(str(epsilon));alpha=mp.acosh(1+e);s=mp.sinh(alpha);total=mp.mpf(0);quiet=0
        for n in range(1,10001):
            k=mp.mpf(2*n+1)
            value=mp.mpf(n*(n+1))/((2*n-1)*(2*n+3))*((2*mp.sinh(k*alpha)+k*mp.sinh(2*alpha))/(4*mp.sinh((mp.mpf(n)+mp.mpf('.5'))*alpha)**2-k*k*s*s)-1)
            total+=value
            quiet=quiet+1 if abs(value)<mp.mpf('1e-32')*max(1,abs(total)) else 0
            if quiet>=8:return float(mp.mpf(4)/3*s*total),n
        raise RuntimeError('Brenner series not converged')

def published_wall_asymptotics(e):
    """Return R normalized using 6*pi*mu*a and q'=(V,aOmega).
    Ytr includes author-fitted linear term; these are near-wall diagnostics,
    not a full finite-gap or independent experimental accuracy certification.
    """
    log=np.log(e)
    return {'Xtt':1/e-.2*log+.9713,
            'Ytt':-8/15*log+.9588,
            'Ytr':4/3*(.1*log+.1895-.4576*e),
            'Xrr':4/3*(1.2021-3*(np.pi*np.pi/6-1)*e),
            'Yrr':4/3*(-.4*log+.3817+1.4578*e),
            'leading_normal':1/e,'GCB_Ytt':-8/15*log+.9588,
            'GCB_Ytr':4/3*(.1*log+.1895),'GCB_Yrr':4/3*(-.4*log+.3817)}
