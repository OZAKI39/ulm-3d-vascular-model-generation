"""P7 distributions, independent RNG streams and absolute cumulative-flux clocks."""
from dataclasses import asdict, dataclass
from copy import deepcopy
import math
import numpy as np
from .sonovue_adapter import read_sonovue
from .rbc_distribution import sample_rbc_geometries
from .rbc import RBCGeometry, GeometryProvenance

H_D = .45
C_MB = 8.5e12
MASTER_SEED = 2026092107
STREAMS = ('MB_SIZE','MB_POSITION','RBC_GEOMETRY','RBC_ORIENTATION','RBC_POSITION','ADMISSION_RETRY')
OLD_RBC_DRIFT = 'DISABLED'


@dataclass(frozen=True)
class ConstantMBConcentrationV0:
    role: str = 'CONSTANT_VALIDATION_WINDOW_V0'
    def number_concentration_m3(self, time_s): return C_MB
    def state(self): return dict(kind='constant', value_m3=C_MB, role=self.role)


class LinearProfile:
    """Continuous piecewise-linear nonnegative Q or C; constant beyond last knot."""
    def __init__(self, times, values):
        self.times = np.asarray(times, float); self.values = np.asarray(values, float)
        if self.times.ndim!=1 or self.values.shape!=self.times.shape or len(self.times)==0 or self.times[0]!=0 or np.any(np.diff(self.times)<=0) or not np.isfinite(self.times).all() or not np.isfinite(self.values).all() or np.any(self.values<0):
            raise ValueError('Finite nonnegative profile with increasing knots starting at zero required')
    def __call__(self, t): return float(np.interp(t, self.times, self.values))
    def number_concentration_m3(self, time_s): return self(time_s)
    def state(self): return dict(kind='linear', times=self.times.tolist(), values=self.values.tolist())


class FluxClock:
    """Exact integral of linear Q times linear C; inversion independent of caller dt."""
    def __init__(self, flow, concentration=None, multiplier=1.):
        self.flow=flow; self.concentration=concentration; self.multiplier=multiplier
        knots = list(flow.times)
        if isinstance(concentration, LinearProfile): knots += list(concentration.times)
        self.knots=np.unique(knots)
        self.coefficients=[]; self.prefix=[0.]
        for i, a in enumerate(self.knots):
            b=self.knots[i+1] if i+1<len(self.knots) else a+1
            q0=flow(a); qs=(flow(b)-q0)/(b-a)
            c0=1. if concentration is None else concentration.number_concentration_m3(a)
            cs=0. if concentration is None else (concentration.number_concentration_m3(b)-c0)/(b-a)
            coef=np.array([q0*c0, q0*cs+qs*c0, qs*cs])*multiplier
            self.coefficients.append(coef)
            if i+1<len(self.knots): self.prefix.append(self.prefix[-1]+self._integral(coef,b-a))
        self.prefix=np.array(self.prefix)
    @staticmethod
    def _integral(coef, dt): return float(dt*(coef[0]+dt*(coef[1]/2+dt*coef[2]/3)))
    def cumulative(self,t):
        if not np.isfinite(t) or t<0: raise ValueError('Nonnegative finite physical time required')
        i=np.searchsorted(self.knots,t,side='right')-1
        return float(self.prefix[i]+self._integral(self.coefficients[i],t-self.knots[i]))
    def time_at(self,target):
        if target<0 or not np.isfinite(target): raise ValueError('Finite nonnegative cumulative threshold required')
        if target==0: return 0.
        first=np.searchsorted(self.prefix,target,side='left')
        if first<len(self.prefix) and self.prefix[first]==target: return float(self.knots[first])
        i=min(np.searchsorted(self.prefix,target,side='right')-1,len(self.knots)-1)
        c=self.coefficients[i]; remainder=target-self.prefix[i]
        if c[1]==c[2]==0:
            return float(self.knots[i]+remainder/c[0]) if c[0]>0 else math.inf
        lo=0.; hi=self.knots[i+1]-self.knots[i]
        for _ in range(64):
            mid=(lo+hi)/2
            if self._integral(c,mid)<remainder: lo=mid
            else: hi=mid
        return float(self.knots[i]+(lo+hi)/2)


def geometry_from_dict(d):
    return RBCGeometry(d['a_m'],d['b_m'],d['c_m'],d['volume_m3'],GeometryProvenance(**d['provenance']))


class PopulationSource:
    def __init__(self, sonovue_root, seed=MASTER_SEED):
        self.sonovue_root=str(sonovue_root); self.seed=seed
        self.rng={name:np.random.default_rng(child) for name,child in zip(STREAMS,np.random.SeedSequence(seed).spawn(len(STREAMS)))}
        self.mb_contract,self.mb_distribution,_=read_sonovue(sonovue_root)
        self.buffer=[]; self.buffer_index=0; self.mb_count=self.rbc_count=0
    def orientation(self, n=1):
        # Uniform S^3 via rotationally invariant Gaussian; antipodal quotient is Haar SO(3).
        q=self.rng['RBC_ORIENTATION'].standard_normal((n,4)); q/=np.linalg.norm(q,axis=1)[:,None]
        return q
    def next_mb(self):
        u=float(self.rng['MB_SIZE'].random()); self.mb_count+=1
        d=float(self.mb_distribution.inverse_cdf(u))
        return dict(species='MB',diameter_um=d,radius_m=d*.5e-6,volume_m3=4*np.pi*(d*.5e-6)**3/3,
            q=[1.,0.,0.,0.],provenance=dict(sequence=self.mb_count,uniform=u,
            sampler_sha256=self.mb_contract['sampler_source_sha256'],histogram_sha256=self.mb_contract['histogram_sha256']))
    def next_rbc(self):
        if self.buffer_index==len(self.buffer):
            seed=int(self.rng['RBC_GEOMETRY'].integers(0,2**63-1))
            pop=sample_rbc_geometries(4096,seed=seed)
            self.buffer=[asdict(RBCGeometry.from_population(pop,i)) for i in range(len(pop.samples))]
            self.buffer_index=0
        g=deepcopy(self.buffer[self.buffer_index]); self.buffer_index+=1; self.rbc_count+=1
        return dict(species='RBC',geometry=g,volume_m3=g['volume_m3'],q=self.orientation()[0].tolist(),
                    provenance=dict(sequence=self.rbc_count,**g['provenance']))
    def state(self):
        return dict(sonovue_root=self.sonovue_root,seed=self.seed,rng={k:deepcopy(v.bit_generator.state) for k,v in self.rng.items()},
            buffer=self.buffer,buffer_index=self.buffer_index,mb_count=self.mb_count,rbc_count=self.rbc_count)
    @classmethod
    def restore(cls,state):
        # No draws and no geometry re-sampling on restore.
        obj=cls.__new__(cls); obj.sonovue_root=state['sonovue_root']; obj.seed=state['seed']; obj.rng={}
        for name,s in state['rng'].items():
            bit=np.random.PCG64(); bit.state=deepcopy(s); obj.rng[name]=np.random.Generator(bit)
        obj.mb_contract,obj.mb_distribution,_=read_sonovue(obj.sonovue_root)
        for k in ('buffer','buffer_index','mb_count','rbc_count'): setattr(obj,k,deepcopy(state[k]))
        return obj


class InjectionScheduler:
    def __init__(self, source, flow, profile=None):
        self.source=source; self.flow=flow; self.profile=profile or ConstantMBConcentrationV0()
        self.mb_clock=FluxClock(flow,self.profile); self.rbc_clock=FluxClock(flow,multiplier=H_D)
        self.time_s=0.; self.next_particle_id=1; self.mb_count=self.rbc_count=0
        self.rbc_volume=0.; self.rbc_sum_compensation=0.; self.next_rbc=source.next_rbc()
    def next_time(self):
        return min(self.mb_clock.time_at(self.mb_count+1),self.rbc_clock.time_at(self.rbc_volume+self.next_rbc['volume_m3']))
    def pop(self):
        tm=self.mb_clock.time_at(self.mb_count+1)
        tr=self.rbc_clock.time_at(self.rbc_volume+self.next_rbc['volume_m3'])
        if tm<=tr:
            event=self.source.next_mb(); self.mb_count+=1; t=tm
        else:
            event=self.next_rbc; self.rbc_count+=1; t=tr
            y=event['volume_m3']-self.rbc_sum_compensation
            total=self.rbc_volume+y; self.rbc_sum_compensation=(total-self.rbc_volume)-y; self.rbc_volume=total
            self.next_rbc=self.source.next_rbc()
        if not np.isfinite(t): raise ValueError('No finite future birth')
        event=deepcopy(event); event.update(particle_id=self.next_particle_id,scheduled_time_s=t,attempt_count=0)
        self.next_particle_id+=1; self.time_s=t
        return event
    def through(self,t):
        if t<self.time_s: raise ValueError('Time cannot decrease')
        while self.next_time()<=t: yield self.pop()
        self.time_s=t
    def state(self):
        return dict(time_s=self.time_s,next_particle_id=self.next_particle_id,mb_count=self.mb_count,rbc_count=self.rbc_count,
            rbc_volume=self.rbc_volume,rbc_sum_compensation=self.rbc_sum_compensation,next_rbc=deepcopy(self.next_rbc),
            mb_cumulative=self.mb_clock.cumulative(self.time_s),rbc_target=self.rbc_clock.cumulative(self.time_s),
            flow=self.flow.state(),profile=self.profile.state(),source=self.source.state())
    @classmethod
    def restore(cls,s):
        obj=cls.__new__(cls); obj.source=PopulationSource.restore(s['source']); obj.flow=LinearProfile(s['flow']['times'],s['flow']['values'])
        p=s['profile']; obj.profile=ConstantMBConcentrationV0() if p['kind']=='constant' else LinearProfile(p['times'],p['values'])
        obj.mb_clock=FluxClock(obj.flow,obj.profile); obj.rbc_clock=FluxClock(obj.flow,multiplier=H_D)
        for k in ('time_s','next_particle_id','mb_count','rbc_count','rbc_volume','rbc_sum_compensation','next_rbc'): setattr(obj,k,deepcopy(s[k]))
        return obj
