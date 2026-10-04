"""Common-strain and component-clamp diagnosis; all outputs are task-local."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import threading
import time

import numpy as np
import psutil
from pd_clot.pfpd_candidate.activation_boundary_control import BoundaryProblem, FULL, AXIAL, solve_intact
from pd_clot.pfpd_candidate import petsc_vi_candidate as vi
from scripts import c7_early_branch as old

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'verification/pfpd_at1_activation_root_cause_and_fix_v1'
RESULTS = ROOT/'results/pfpd_at1_activation_root_cause_and_fix_v1'
STRAINS = (.0025, .003)


def save(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False,
                               default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v.item())+'\n')


def read(path):
    return json.loads(Path(path).read_text())


def arrays(path):
    with np.load(path, allow_pickle=False) as z:
        return {k:z[k] for k in z.files}


def declare():
    OUT.mkdir(parents=True, exist_ok=True); RESULTS.mkdir(parents=True, exist_ok=True)
    path = OUT/'CONFIG.json'
    if path.exists():
        return read(path)
    mem = psutil.virtual_memory()
    workers = max(1, min(6, psutil.cpu_count(logical=False)-2,
                        int((mem.available-.15*mem.total-.75*1024**3)/(.65*1024**3))))
    p = BoundaryProblem(5,0)
    config = dict(task='PFPD_AT1_ACTIVATION_ROOT_CAUSE_AND_FIX_V1',
        model='PFPD_AT1_SPECTRAL_CANDIDATE_V1', paper_exact_equivalence_is_gate=False,
        cleanup_byte_failures_are_gate=False, phase_A_strains=STRAINS, levels=[5,7],
        physics=p.config, L_m=old.old.v2.L, delta_m=p.family.horizon, Yc_Pa=p.yc,
        material=dict(mu=p.material.mu, bulk=p.material.bulk, lame=p.material.lame),
        workers=workers, physical_cores=psutil.cpu_count(logical=False), logical_cores=psutil.cpu_count(),
        RAM_total=mem.total, RAM_available_at_start=mem.available, minimum_headroom_fraction=.15,
        MPI_ranks=1, BLAS_threads=1, GPU='GPU_NOT_USED_BY_DESIGN',
        petsc_settings=old.prior.config()['primary'], force_tolerance=1e-9,
        reaction_balance_tolerance=1e-8, positive_J_required=True, boundary_exact_required=True,
        affine_oracle='F=diag(1,1,1+epsilon), with bottom axial translation origin; nu=0. Diagnostic, not presumed equilibrated.',
        axial_selection='abs(n_z)>=0.95', first_layer='axial bond midpoint d_clamp/h=0.5, absolute tolerance 1e-8',
        A_psi_denominator_floor_fraction_Yc=1e-12,
        affine_uniform_CV_tolerance=1e-8, affine_uniform_relative_span_tolerance=1e-8,
        equilibrium_hotspot_Apsi_threshold=1.2,
        profile=dict(primary='arithmetic mean A_psi of axial bonds at each distinct midpoint distance',
                     supplemental=['Psi_plus/Yc','q/Yc','P95 A_psi','max A_psi'],
                     interpolation='1001 fixed equally spaced samples on overlap of the two profile domains; linear, no extrapolation',
                     common_grid_points=1001, strong_ratio=.5, require_both_strains=True),
        clamp=dict(hotspot_metric='maximum first-layer axial Psi_plus/Yc',
                   discrepancy_metric='abs(max_driver_C7-max_driver_C5)/mean(max_driver_C7,max_driver_C5)',
                   strong_hotspot_reduction=.20, strong_discrepancy_reduction=.30, weak_change=.05),
        support=dict(metrics=['mean endpoint omega0','max endpoint shape condition'],
                     strong_abs_spearman=.7, weak_abs_spearman=.3,
                     definition='STRONG if same support metric has |rho|>=.7 at both levels; WEAK if all |rho|<.3; otherwise MIXED. Association only, not causation.',
                     bulk='most-interior available points, never a full-horizon bulk'),
        activation_bracket_width=1e-7, activation_discrepancy_denominator='mean of C5,C7',
        candidate_selection='Only after Phase A decision. One physical change; original material, horizon, phase law and PETSc settings retained.',
        early_phase_default='NOT_RUN; optional, only after clear candidate improvement',
        historical_writers_invoked=False, new_repository_manifests=False)
    save(path,config)
    return config


def resource_monitor(stop, started):
    proc=psutil.Process(); peak=0; minimum=psutil.virtual_memory().available; samples=0
    while not stop.is_set():
        try:
            members=[proc]+proc.children(recursive=True)
            rss=sum(q.memory_info().rss for q in members if q.is_running())
            available=psutil.virtual_memory().available
            peak=max(peak,rss);minimum=min(minimum,available);samples+=1
            save(OUT/'RESOURCES.json',dict(elapsed_seconds=time.monotonic()-started,tree_peak_RSS_bytes=peak,
                minimum_observed_available_bytes=minimum,minimum_observed_headroom_fraction=minimum/psutil.virtual_memory().total,
                samples=samples,workers=read(OUT/'CONFIG.json')['workers'],GPU='GPU_NOT_USED_BY_DESIGN'))
        except psutil.NoSuchProcess:
            pass
        stop.wait(2.)


def gates_and_summary(p,x,s):
    state=p.state(x,s); m=p.metrics(x,s); Y=state.W0; q=2*(p.yc-Y)
    J=np.linalg.det(state.Ftilde);Jbar=np.linalg.det(state.Fbar)
    force=state.force
    top=-force[p.top].sum(axis=0);bottom=-force[p.bottom].sum(axis=0)
    scale=max(np.linalg.norm(top),np.linalg.norm(bottom),p.config['mechanical_equilibrium']['force_floor_N'])
    reaction=np.where(p.fixed_dof,-force,0.)
    reaction_balance=float(np.linalg.norm(reaction.sum(axis=0))/scale)
    safety=(np.isfinite(x).all() and np.linalg.norm(x-p.boundary,axis=1).max()<=old.old.GUARD['safety']['maximum_displacement_m']
            and np.linalg.norm(x[:,1:3],axis=1).max()<=old.old.GUARD['pipe']['radius_m']*old.old.GUARD['safety']['maximum_radius_ratio'])
    gates=dict(force=m['force_residual']<=1e-9,positive_J=bool(np.all(J>0) and np.all(Jbar>0)),
               finite_J=bool(np.isfinite(J).all() and np.isfinite(Jbar).all()),safety=bool(safety),
               BC_exact=bool(np.array_equal(x[p.fixed_dof],p.boundary[p.fixed_dof])),
               phase_exact_zero=bool(np.count_nonzero(s)==0),reaction_balance=reaction_balance<=1e-8,
               rigid_modes_removed=p.gauge_rank()==6)
    idx=int(np.argmin(q));mid=p.family.X[p.family.pairs].mean(axis=1)
    d=old.old.v2.L/2-abs(mid[:,2])
    row=dict(level=p.n,strain=p.strain,BC_MODE=p.mode,metrics=m,minimum_J=float(J.min()),
             minimum_Jbar=float(Jbar.min()),max_Psi_plus_Pa=float(Y.max()),P95_Psi_plus_Pa=float(np.percentile(Y,95)),
             min_q_Pa=float(q.min()),max_driver_normalized=float(Y.max()/p.yc),
             reaction_top_N=top,reaction_bottom_N=bottom,reaction_balance_relative=reaction_balance,
             reaction_stress_Pa=float(top[2]/old.old.v2.L**2),energy_J=state.energy,
             min_q_bond_ID=idx,min_q_midpoint_m=mid[idx],min_q_d_clamp_m=float(d[idx]),
             min_q_d_over_h=float(d[idx]/p.cloud.spacing),min_q_d_over_delta=float(d[idx]/p.family.horizon),
             gauge_constraints=p.gauge,free_components=p.nx,gates=gates,accepted=all(gates.values()))
    return row


def equilibrium(n,e,mode,pool,label=None,initial=None):
    label=label or f'C{n}/{mode}/{e:.8f}'
    folder=OUT/'equilibria'/label;p=BoundaryProblem(n,e,mode)
    if (folder/'SUMMARY.json').exists():
        row=read(folder/'SUMMARY.json');a=arrays(folder/'STATE.npz')
        assert row['level']==n and row['strain']==e and row['BC_MODE']==mode
        assert gates_and_summary(p,a['x'],a['s'])['accepted']
        return row,a
    source=None
    if n==7 and mode==FULL and e in STRAINS:
        source=old.OUT/f'intact/scan_{1 if e==.0025 else 2:02d}'
    if source:
        oldrow=read(source/'SUMMARY.json');assert oldrow['strain']==e and oldrow['intact_accepted']
        a=arrays(source/'STATE.npz');x=a['x'];s=a['s']
        solver=dict(method='REUSED_EXACT_EXISTING_EQUILIBRIUM',source=str(source.relative_to(ROOT)),
                    source_state_sha256=hashlib.sha256((source/'STATE.npz').read_bytes()).hexdigest(),new_solve=False)
    else:
        if initial is None:
            x=p.affine()
        else:
            x=initial['x'].copy()
            x[:,2]+=(e-float(initial['strain']))*(p.family.X[:,2]+old.old.v2.L/2)
        assert psutil.virtual_memory().available/psutil.virtual_memory().total>=.15,'RAM_HEADROOM'
        x,s,solver=solve_intact(p,x,pool,read(OUT/'CONFIG.json')['workers'])
        solver['new_solve']=True
    row=gates_and_summary(p,x,s);row['solver']=solver
    save(folder/'SUMMARY.json',row)
    if row['accepted']:
        np.savez_compressed(folder/'STATE.npz',x=x,s=s,strain=e)
    else:
        np.savez_compressed(folder/'FAILED_STATE.npz',x=x,s=s,strain=e)
    print('ACCEPTANCE',n,mode,e,row['accepted'],row['metrics']['force_residual'],row['max_driver_normalized'],flush=True)
    assert row['accepted'],'MECHANICAL_GATES_FAILED: '+label
    return row,dict(x=x,s=s,strain=e)


def phase_a(pool):
    rows=[]
    for mode in (FULL,AXIAL):
        for n in (5,7):
            prev=None
            for e in STRAINS:
                row,prev=equilibrium(n,e,mode,pool,initial=prev)
                rows.append(row);save(OUT/'COMMON_STRAIN_EQUILIBRIA.json',rows)
    return rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['phase_a'])
    args=parser.parse_args();config=declare();stop=threading.Event();started=time.monotonic()
    monitor=threading.Thread(target=resource_monitor,args=(stop,started),daemon=True);monitor.start()
    pool=vi.make_pool(config['workers'])
    try:
        if args.stage=='phase_a':phase_a(pool)
    finally:
        if pool:pool.shutdown()
        stop.set();monitor.join()


if __name__=='__main__':
    main()
