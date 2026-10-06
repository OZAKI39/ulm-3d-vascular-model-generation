"""Versioned R2 motion: S0 diagnostic and converged S1 semi-implicit Euler.

S1 freezes F_acoustic and u_flow at the beginning of each substep. It solves
 A(x0+dt U) U = F0 + D u_flow0 + J(x0+dt U)' lambda.
Contact inequalities are linearized about the current endpoint and expressed
in the ONE x0 displacement anchor. Historical tangent rows are not stacked.
The resistance is never differentiated as a position energy.
"""
from collections import OrderedDict
from time import perf_counter
import numpy as np
from scipy import sparse
from scipy.spatial.distance import pdist
from .constrained_motion import (ConstrainedLubrication, Constraints, sparse_rows,
                                  MotionSolverFailure, SubdivideInterval, kkt_residual)
from .translation_backends import _pairs, _gap_vectors
from .motion_qp_r2 import NativeQP

BACKEND='constrained_lubrication_r2'


class MotionR2(ConstrainedLubrication):
    name=BACKEND

    def __init__(self, viscosity_Pa_s, *, solver='prox', **kwargs):
        super().__init__(viscosity_Pa_s, **kwargs)
        self.qp_factory=lambda: NativeQP(solver)
        self.qp=self.qp_factory();self.qp_solver=solver
        self.method='STAGE5A_R2_S1_SEMI_IMPLICIT_EULER_V1'
        self.stage_workspaces=OrderedDict()

    def constrained_stage(self,result,constraints,ids,stage,feasible=None):
        if stage not in self.stage_workspaces:
            self.stage_workspaces[stage]=OrderedDict()
        self.stage_workspaces.move_to_end(stage)
        self.component_workspaces=self.stage_workspaces[stage]
        answer=super().constrained(result,constraints,ids,feasible=feasible)
        while len(self.component_workspaces)>128:
            self.component_workspaces.popitem(last=False)
        while len(self.stage_workspaces)>12:
            self.stage_workspaces.popitem(last=False)
        return answer


def anchored_constraints(x0,xk,r,wall0,wallk,dt,ids,tolerance,reach_velocity):
    """Every pair potentially reached by the known candidate plus activation band.

    Broad phase is augmented by both endpoints; final all-pair path check rejects
    any omission. Initial tolerance-band gaps are retained, never clamped in data.
    """
    radius=r[0]
    reach=2*np.max(np.linalg.norm(reach_velocity-reach_velocity.mean(0),axis=1),initial=0)*dt
    pairs=np.unique(np.vstack([_pairs(x0,2.1*radius+reach),_pairs(xk,2.1*radius)]),axis=0)
    hk,nk=_gap_vectors(xk,pairs,radius)
    h0,_=_gap_vectors(x0,pairs,radius)
    wg0,_,valid0=wall0;wg,nw,valid=wallk
    if not valid.all() or not valid0.all():
        raise MotionSolverFailure(dict(reason='GEOMETRY_DATA_CORRUPTED'))
    wi=np.flatnonzero(np.minimum(wg,wg0)<=np.maximum(.1*radius,np.linalg.norm(reach_velocity,axis=1)*dt))
    g=sparse_rows(pairs,nk,wi,nw[wi],len(r))
    initial=np.r_[h0,wg0[wi]];current=np.r_[hk,wg[wi]]
    allowed=np.where(initial<=tolerance,initial,0.)
    lower=(allowed-current+g@(xk-x0).ravel())/dt
    keys=tuple(('pair',int(ids[i]),int(ids[j])) for i,j in pairs)
    keys+=tuple(('wall',int(ids[i])) for i in wi)
    return Constraints(g,lower,keys,pairs)


class SemiImplicit:
    def __init__(self,backend,geometry,forcing,sweep,*,method='S1',tolerance=1e-11,
                 max_iterations=32,nonlinear_rtol=1e-9,atol_R=.002,rtol=.01,gap_rtol=.02,
                 relaxation=1.):
        self.backend,self.geometry,self.forcing,self.sweep=backend,geometry,forcing,sweep
        self.method,self.tolerance=method,tolerance
        self.max_iterations,self.nonlinear_rtol=max_iterations,nonlinear_rtol
        self.atol_R,self.rtol,self.gap_rtol=atol_R,rtol,gap_rtol
        self.relaxation=relaxation
        self.timings={k:0. for k in ['forcing','wall','assembly','qp','geometry','error']}
        self.calls=0

    def wall(self,x,r):
        start=perf_counter();w=self.geometry.wall_gap_and_normal(x,r)
        self.timings['wall']+=perf_counter()-start
        return w

    def advance(self,x,r,ids,t,dt,stage):
        start=perf_counter();f=self.forcing(x,r,t)
        self.timings['forcing']+=perf_counter()-start
        force=f['primary_force_N']+f['secondary_force_N'];flow=f['flow_velocity_m_s']
        w0=self.wall(x,r);xk=x.copy();uk=np.zeros_like(x)
        diagnostics=[];last_delta=None
        for iteration in range(1 if self.method=='S0' else self.max_iterations):
            wk=self.wall(xk,r)
            start=perf_counter();a=self.backend.solve(xk,r,force,flow,wk)
            self.timings['assembly']+=perf_counter()-start
            c=anchored_constraints(x,xk,r,w0,wk,dt,ids,self.tolerance,a.velocity)
            start=perf_counter()
            u,cache,diag=self.backend.constrained_stage(a,c,ids,stage,feasible=uk.ravel())
            self.timings['qp']+=perf_counter()-start;self.calls+=1
            y=x+dt*u
            delta=np.max(np.linalg.norm(u-uk,axis=1),initial=0)
            scale=max(np.max(np.linalg.norm(u,axis=1),initial=0),1e-12)
            diag.update(linear_residual=a.diagnostics['relative_residual'], nonlinear_delta=delta/scale)
            diagnostics.append(diag)
            if self.method=='S0' or delta/scale<=self.nonlinear_rtol:
                # Re-evaluate the unlagged endpoint equations with the returned U,
                # multipliers and the SAME fixed Stage5A normalization.
                if self.method=='S1':
                    we=self.wall(y,r);beg=perf_counter();ae=self.backend.solve(y,r,force,flow,we)
                    self.timings['assembly']+=perf_counter()-beg
                    ce=anchored_constraints(x,y,r,w0,we,dt,ids,self.tolerance,u)
                    vs=max(float(np.max(abs(ae.velocity))),1e-30)
                    an=ae.matrix/ae.scale;bn=an@(ae.velocity.ravel()/vs)
                    lam=np.array([cache.dual.get(key,0.)/(ae.scale*vs) for key in ce.ids])
                    k,terms=kkt_residual(an,bn,ce.matrix,ce.lower/vs,u.ravel()/vs,lam)
                    diag['nonlinear_kkt']=k;diag['nonlinear_kkt_terms']=terms
                    if k>1e-9:
                        uk=u;xk=y
                        continue
                    diag['endpoint_linear_residual']=ae.diagnostics['relative_residual']
                else:
                    diag['nonlinear_kkt']=None
                begin=perf_counter();safety=self.sweep(x,y,r,self.tolerance)
                self.timings['geometry']+=perf_counter()-begin
                if not safety['pass_gate']:
                    raise SubdivideInterval('SWEEP_NOT_CERTIFIED',dict(dt_s=dt,safety=safety))
                return y,u,cache,dict(qp=diag,iterations=iteration+1,
                    kkt=max(z['kkt'] for z in diagnostics),
                    linear_residual=max(z['linear_residual'] for z in diagnostics),
                    qp_iterations=sum(z['iterations'] for z in diagnostics),
                    fallbacks=sum(s.get('fallbacks',0) for z in diagnostics for s in z['solves']),safety=safety)
            # Relaxation is ONLY nonlinear iteration, not a force or motion term.
            alpha=self.relaxation
            if last_delta is not None and delta>last_delta*1.2:
                alpha=min(alpha,.5)
            uk=(1-alpha)*uk+alpha*u;xk=x+dt*uk;last_delta=delta
        raise SubdivideInterval('NONLINEAR_NOT_CONVERGED',dict(dt_s=dt,iterations=len(diagnostics),
                                    last_delta=diagnostics[-1]['nonlinear_delta']))

    def trial(self,x,r,ids,t,dt):
        full,_,_,df=self.advance(x,r,ids,t,dt,'full')
        half,uh,ch,dh=self.advance(x,r,ids,t,dt/2,'half1')
        fine,u,c,d=self.advance(half,r,ids,t+dt/2,dt/2,'half2')
        begin=perf_counter()
        diff=fine-full
        position=float(np.max(np.linalg.norm(diff,axis=1),initial=0))
        relative=float(np.max(np.linalg.norm(diff-diff.mean(0),axis=1),initial=0))
        displacement=np.linalg.norm(fine-x,axis=1).max(initial=0)
        posbudget=self.atol_R*r[0]+self.rtol*displacement
        h0=pdist(x)-2*r[0];hf=pdist(full)-2*r[0];hh=pdist(fine)-2*r[0]
        pairratio=abs(hh-hf)/(self.atol_R*r[0]+self.gap_rtol*np.maximum(h0,.001*r[0]))
        w0=self.wall(x,r)[0];wf=self.wall(full,r)[0];wh=self.wall(fine,r)[0];valid=np.isfinite(w0)
        wallratio=abs(wf[valid]-wh[valid])/(self.atol_R*r[0]+self.gap_rtol*np.maximum(w0[valid],.001*r[0]))
        ratios=dict(position=position/posbudget,relative=relative/posbudget,
                    pair=float(pairratio.max(initial=0)),wall=float(wallratio.max(initial=0)))
        error=max(ratios.values());self.timings['error']+=perf_counter()-begin
        diag=dict(error_ratio=error,error_terms=ratios,position_error_m=position,
                  kkt=max(q['kkt'] for q in [df,dh,d]),
                  linear_residual=max(q['linear_residual'] for q in [df,dh,d]),
                  nonlinear_iterations=sum(q['iterations'] for q in [df,dh,d]),
                  qp_iterations=sum(q['qp_iterations'] for q in [df,dh,d]),
                  fallbacks=sum(q['fallbacks'] for q in [df,dh,d]))
        if error>1:
            raise SubdivideInterval('LOCAL_ERROR',dict(diag,dt_s=dt))
        # Actual adopted path consists of TWO certified segments, both recorded.
        return [(half,uh,ch,dh),(fine,u,c,d)],diag
