"""Explicit opt-in low-pressure intervals; historical arrays/queries preserved."""
import numpy as np
from scripts.straight_tube_exp1.f10ad_multifrequency.interface import FixedResponse
from .common import BASE,sha
class ExtendedResponse:
    def __init__(self,old,physics):
        low=FixedResponse(BASE/'low_pressure/diagnostic_response.h5',physics)
        self.old=old;self.low=low;boundary=float(old.knots[0]);ix=low.knots<boundary
        extra=low.knots[ix];vv=low.knot_values[ix]
        if len(extra)<1:raise ValueError('No validated low-pressure intervals')
        # Preserve exactly the historical endpoint value and all historical intervals.
        joined=np.concatenate([vv,old.knot_values[:1]])
        self.knots=np.r_[extra,old.knots]
        self.knot_values=np.concatenate([vv,old.knot_values])
        self.knot_valid=np.r_[np.ones(len(extra),np.int32),old.knot_valid]
        self.values=np.concatenate([np.stack([joined[:-1],joined[1:]],1),old.values])
        self.valid=np.r_[np.ones(len(extra),np.int32),old.valid]
        self.lower=self.knots[:-1];self.upper=self.knots[1:]
        self.metadata=dict(model='F10AE_VALIDATED_LOW_PRESSURE_EXTENSION',historical=old.metadata,low=low.metadata,low_sha256=sha(BASE/'low_pressure/diagnostic_response.h5'),interval_min_Pa=float(self.knots[0]),interval_max_Pa=float(self.knots[-1]),historical_values_preserved=True,extrapolation=False)
    def query(self,p):
        p=np.asarray(p);out=np.zeros(p.shape+(8,),complex);a=abs(p)
        if not np.isfinite(p).all() or np.any((a>0)&(a<self.knots[0])):raise ValueError('Outside validated low-pressure range; no clipping')
        high=a>=self.old.knots[0];low=(a>0)&~high
        if np.any(high):out[high]=self.old.query(p[high])
        if np.any(low):
            n=int(np.count_nonzero(self.knots<self.old.knots[0]));kn=self.knots[:n+1];v=self.knot_values[:n+1]
            for k in range(8):out[...,k][low]=(np.interp(a[low],kn,v[:,k].real)+1j*np.interp(a[low],kn,v[:,k].imag))*np.exp(1j*(k+1)*np.angle(p[low]))
        return out
