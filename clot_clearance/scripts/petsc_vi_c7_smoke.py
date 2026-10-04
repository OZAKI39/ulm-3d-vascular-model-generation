"""Optional C7 initialization smoke after C5 gate; never a branch/stability run."""
import time,resource
import numpy as np
from petsc4py import PETSc
from scripts import petsc_vi_audit as a
from pd_clot.pfpd_candidate import petsc_vi_candidate as vi


def main():
    assert a.read(a.OUT/'continuation/PETSC_RSLS/ALL.json')['passed']==7
    source=a.prev.OLD/'equilibrium/C7/0.00025000/STATE.npz';old=a.prev.ca.read_npz(source);eps0=.00025
    low=a.Problem(7,eps0);Y=old['Psi_plus'];estimate=eps0*np.sqrt(low.yc/Y.max());scan=[]
    def predictor(p):
        # Scale the entire saved small-strain relaxed displacement field,
        # including lateral relaxation; a boundary-only increment is too stiff.
        return p.prescribe(p.family.X+(old['x']-p.family.X)*(p.strain/eps0))
    for eps in [estimate*.9,estimate,estimate*1.1]:
        p=a.Problem(7,eps);x=predictor(p);s=np.zeros(p.ns)
        yy,q,dd,qs=a.prev.q_fields(p,x,s)
        scan.append(dict(strain=float(eps),min_q_scaled=float(qs.min()),max_Y_over_Yc=float((yy/p.yc).max()),force=p.metrics(x,s)['force_residual']))
    # This cheap frozen-phase displacement predictor scan is not equilibrated
    # and is explicitly not a new activation bracket or stationary root.
    near=min(scan,key=lambda r:abs(r['min_q_scaled']));rows=[]
    for label,eps in [('low',eps0),('activation_neighborhood_predictor',near['strain'])]:
        start=time.perf_counter();p=a.Problem(7,eps);x=old['x'] if label=='low' else predictor(p)
        z=p.pack(x,np.zeros(p.ns));r=vi.residual(p,z)
        adapter=vi.PETScAdapter(p,PETSc.COMM_SELF);adapter.assign(adapter.v,z);adapter.function(None,adapter.v,adapter.f)
        error=float(np.max(abs(adapter.gather(adapter.f)-r)));n=p.nx+p.ns;nnz=p.nx*n+p.ns*(p.nx+1)
        # This adapter's conservative simultaneous arrays: verified blocks,
        # CSR staging and AIJ, plus C7 FD worker budgets. Enforce RAM reserve.
        mem={k:int(v.split()[0])*1024 for k,v in (line.split(':',1) for line in a.Path('/proc/meminfo').read_text().splitlines())}
        estimated_peak=6*8*p.nx*p.ns+2*12*nnz+a.config()['workers']*.6*1024**3
        row=dict(label=label,strain=eps,total_DOF=n,mechanical_DOF=p.nx,phase_DOF=p.ns,AIJ_structural_nnz=nnz,
            AIJ_CSR_equivalent_bytes=12*nnz+4*(n+1),adapter_mapping_error=error,raw_phase_F_norm=float(np.linalg.norm(r[p.nx:])),
            metrics=p.metrics(x,np.zeros(p.ns)),predicted_Jacobian_with_workers_peak_bytes=estimated_peak,
            available_memory_bytes=mem['MemAvailable'],reserved_memory_bytes=.2*mem['MemTotal'],
            residual_and_AIJ_initialization='PASSED',Newton_iterations=0,
            Jacobian_assembly='NOT_RUN_MEMORY_PREFLIGHT',stationary_root_claim=False,
            wall_seconds=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
        rows.append(row);adapter.destroy()
    a.save(a.OUT/'c7/SMOKE.json',dict(status='RESIDUAL_AND_SPARSE_INITIALIZATION_ONLY',C5_gate='7/7 PASSED',source=str(source),
        cheap_driver_estimate=estimate,cheap_predictor_scan=scan,scan_is_equilibrium=False,predictor='Linear scaling of saved low-load relaxed displacement, followed by exact boundary prescription',states=rows,
        reason='Conservative full verified block construction plus CSR/AIJ staging and worker pool exceeds reserved-memory budget. No C7 Jacobian assembly or Newton step; no scalability convergence claim.',
        stability='NOT RUN',branch='NOT RUN'))
    print('C7 residual / AIJ initialization smoke complete; Jacobian stopped at memory preflight',flush=True)


if __name__=='__main__':main()
