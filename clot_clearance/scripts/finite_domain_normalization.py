"""Algebraic localization after a failed one-horizon bulk test; no BC/model change."""
import resource,time
import numpy as np
from scripts.finite_domain_boundary_audit import CASES,OUT,RESULTS,DELTA,EPS,material,save,table
from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.boundary_force_diagnostics import (
    make_cube,continuum,force_terms,endpoint_weight,self_neighbor_terms,full_lattice_omega0)


def main():
    start=time.monotonic();summary=[]
    for name,n,L in CASES:
        c=make_cube(n,L,DELTA);f=c.family;m=material();ref=continuum(m,L,EPS)
        x=f.X@ref['F'].T;x[:,2]+=EPS*L/2;st=evaluate(f,x,m);core=c.full_support
        channels=self_neighbor_terms(f,st)
        norm=np.linalg.norm(st.force)
        assert np.linalg.norm(sum(channels.values())-st.force)/norm<1e-12
        omega=full_lattice_omega0(c.h,DELTA);wref=f.omega/omega
        uniform=force_terms(f,st,wref)['total'];excess=force_terms(f,st,f.wbar-wref)['total']
        trunc=force_terms(f,st,endpoint_weight(f,~core))['total'];full=force_terms(f,st,endpoint_weight(f,core))['total']
        assert np.linalg.norm(uniform+excess-st.force)/norm<1e-12
        assert np.linalg.norm(trunc+full-st.force)/norm<1e-12
        rows=[]
        first=(c.ijk[:,0]==n//2)&(c.ijk[:,1]==n//2)&np.isin(c.ijk[:,2],[1,n-2])
        for pid in np.flatnonzero(core|first):
            d=float(st.force[pid,2]);active=abs(d)>1e-12*ref['reaction_N']
            row=dict(case=name,ID=int(pid),Z_m=f.X[pid,2],d_min_over_delta=c.distances[pid].min()/DELTA,
                     full_support=bool(core[pid]),first_interior_center=bool(first[pid]),force_z_N=d)
            for k,v in {**channels,'uniform_wbar_component':uniform,'normalization_excess_component':excess,
                        'truncated_endpoint_component':trunc,'full_support_endpoint_component':full}.items():
                row[k+'_z_N']=float(v[pid,2]);row[k+'_signed_fraction']=float(v[pid,2]/d) if active else None
            rows.append(row)
        table(RESULTS/'tables'/f'{name}_coefficient_channels.csv',rows)
        summary.append(dict(case=name,lattice_omega0_m3=omega,actual_omega0_min_ratio=float(f.omega0.min()/omega),actual_omega0_max_ratio=float(f.omega0.max()/omega),
            channel_reconstruction_relative_error=float(np.linalg.norm(sum(channels.values())-st.force)/norm),
            weight_reconstruction_relative_error=float(np.linalg.norm(uniform+excess-st.force)/norm),
            truncated_split_relative_error=float(np.linalg.norm(trunc+full-st.force)/norm),
            meaning='Exact algebraic summands, never an alternate constitutive or boundary model',rows=rows))
    save(OUT/'NORMALIZATION_CHANNELS.json',dict(cases=summary,wall_seconds=time.monotonic()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))


if __name__=='__main__':main()
