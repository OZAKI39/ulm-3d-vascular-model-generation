"""Unmodified PyStokes wallBounded API, direct SI six-impulse extraction."""
import numpy as np
import pystokes

def get_wall_mobility(radius_m, gap_m, viscosity_pa_s):
    a,h,eta=map(float,(radius_m,gap_m,viscosity_pa_s))
    if min(a,h,eta)<=0: raise ValueError('Strictly positive radius, surface gap and viscosity required')
    rbm=pystokes.wallBounded.Rbm(radius=a,particles=1,viscosity=eta)
    r=np.array([0.,0.,a+h]);M=np.zeros((6,6))
    for k in range(6):
        F=np.zeros(3);T=np.zeros(3);v=np.zeros(3);o=np.zeros(3)
        (F if k<3 else T)[k%3]=1.
        rbm.mobilityTT(v,r,F);rbm.mobilityTR(v,r,T)
        rbm.mobilityRT(o,r,F);rbm.mobilityRR(o,r,T)
        M[:,k]=np.r_[v,o]
    return M,{'implementation':'pystokes.wallBounded.Rbm mobilityTT/TR/RT/RR','version':pystokes.__version__,'native_units':'direct SI','native_order':'N=1: x,y,z for every vector; general N structure-of-arrays','six_impulses':'1 N in F columns; 1 N*m in T columns','output_accumulation':'fresh zero v/o per column','sign_correction_applied':False,'matrix_symmetrized':False,'all_blocks_available':True}
