"""Frozen nonaffine localization audit; independent CPU oracle jobs, serial solves."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gc
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import resource
import threading
import time
import numpy as np
import psutil
from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.energy import Material
from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import block_family
from pd_clot.pfpd_candidate.nonlocal_boundary.diagnostics import support_profile
from pd_clot.pfpd_candidate.nonlocal_boundary.mechanics_adapter import solve
from pd_clot.pfpd_candidate.nonlocal_boundary.nonaffine_manufactured import QuadraticAxial
from pd_clot.pfpd_candidate.nonlocal_boundary.bulk_reference import ManufacturedStencil, oracle_job
from pd_clot.pfpd_candidate.nonlocal_boundary.localized_conversion import localize
from pd_clot.pfpd_candidate.nonlocal_boundary.residual_decomposition import balances

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'verification/nonlocal_boundary_nonaffine_localization_v1'
OLD = ROOT/'verification/nonlocal_boundary_formulation_v1'


def plain(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    raise TypeError(type(x))


def save(name, data): (OUT/name).write_text(json.dumps(data, indent=2, default=plain, allow_nan=False)+'\n')
def load(name): return json.loads((OUT/name).read_text())
def csvsave(name, rows):
    with (OUT/name).open('w', newline='') as s:
        w = csv.DictWriter(s, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def intact_material(): return Material(model='svk', mu=1000., bulk=666.6666666666666)
def normmax(a): return float(np.linalg.norm(a, axis=-1).max())
def integrity():
    hashes = load('PROTECTED_INPUT_HASHES.json')
    different = [p for p, h in hashes.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != h]
    return dict(protected_files=len(hashes), differences=different, PASS=not different)


def require(value, label):
    if not value:
        save('HARD_STOP.json', dict(reason=label, no_threshold_or_physics_changes=True))
        raise RuntimeError('HARD STOP: '+label)


def oracle_stage(cfg):
    field = QuadraticAxial(); h = cfg['h_m']; L = cfg['L_m']; scale = field.scales(h); Fref = scale['F_ref_N']
    continuum = field.continuum(); r = field.response(np.linspace(-L/2, L/2, 19))
    exact = dict(field=dict(epsilon=field.epsilon, beta=field.beta, L_m=L), scales=scale,
                 Z_m=np.linspace(-L/2, L/2, 19), response=r, continuum=continuum,
                 provenance='Independent formulas in nonaffine_manufactured.py; local SVK lambda=0',
                 VERIFIED=bool(continuum['force_closure_scaled']<=cfg['gate'] and continuum['max_work_closure_scaled']<=cfg['gate']))
    save('EXACT_NONAFFINE_REFERENCE.json', exact); save('NONAFFINE_AXIAL_EXACT_REFERENCE.json', exact)
    require(exact['VERIFIED'], 'local continuum reference')
    # Evaluate exact D1/D2/D3 graph closure once per distinct scaffold geometry.
    closure = []
    for n, offset in [(27, [0, 0, 0]), (29, [1, -1, 1])]:
        s = ManufacturedStencil(field, n); c = s.closure(s.id_at_offset(offset)); c.update(n=n, offset=offset)
        closure.append(c)
    save('BULK_STENCIL_CLOSURE.json', closure)
    require(all(c['FULL_FORCE_STENCIL_CLOSED'] for c in closure), 'bulk dependency graph/volumes')
    del s; gc.collect()
    zs = np.linspace(-L/2, L/2, cfg['N'])
    specs = [dict(label='primary', level=k, Z=float(z)) for k, z in enumerate(zs)]
    for k in [0, 9, 18]:
        specs += [dict(label='affine', level=k, Z=float(zs[k]), beta=0.),
                  dict(label='shifted', level=k, Z=float(zs[k]), n=29, offset=[1, -1, 1], x=3*h, y=-2*h),
                  dict(label='reflected', level=k, Z=float(-zs[k]), beta=-field.beta)]
    rows = []
    # Spawn avoids inheriting an expensive dependency-closure scratch heap.
    with ProcessPoolExecutor(cfg['oracle_workers'], mp_context=mp.get_context('spawn')) as pool:
        for row in pool.map(oracle_job, specs):
            rows.append(row); print('ORACLE', row['label'], row['level'], row['q_bulk_N_m3'], flush=True)
    save('BULK_ORACLE_RAW.json', rows)
    primary = [r for r in rows if r['label']=='primary']; q = np.array([r['q_bulk_N_m3'] for r in primary])
    b = field.response(zs)['b']; error = q+b
    checks = []
    for row in rows[19:]:
        k = row['level']; actual = np.array(row['q_bulk_N_m3'])
        target = np.zeros(3) if row['label']=='affine' else q[k].copy()
        if row['label']=='reflected': target[2] *= -1
        err = float(np.linalg.norm(actual-target)*h**3/Fref)
        checks.append(dict(label=row['label'], level=k, local_force_scaled_error=err, PASS=err<=cfg['gate']))
    # All transverse components should vanish for a translationally uniform axial field.
    transverse = float(np.linalg.norm(q[:, :2], axis=1).max()*h**3/Fref)
    verified = all(c['PASS'] for c in checks) and transverse <= cfg['gate']
    summary = dict(checks=checks, transverse_force_scaled=transverse, VERIFIED=verified,
                   max_local_body_N_m3=normmax(b), max_bulk_q_N_m3=normmax(q), max_model_error_N_m3=normmax(error),
                   model_error_relative_to_body=normmax(error)/scale['b_ref_N_m3'],
                   model_error_relative_to_stress_gradient=normmax(error)/scale['stress_gradient_scale_N_m3'])
    summary['BULK_MODEL_DIFFERENCE'] = ('SMALL' if summary['model_error_relative_to_body']<=cfg['bulk_difference_threshold'] else 'MATERIAL') if verified else 'INCONCLUSIVE'
    save('BULK_ORACLE_VALIDATION.json', summary)
    csvsave('BULK_BAQP_REFERENCE.csv', [dict(level=k, Z_m=float(z), **{f'q_bulk_{c}_N_m3':q[k,j] for j,c in enumerate('xyz')},
                                               **{f'b_local_{c}_N_m3':b[k,j] for j,c in enumerate('xyz')},
                                               **{f'e_model_{c}_N_m3':error[k,j] for j,c in enumerate('xyz')}) for k,z in enumerate(zs)])
    require(verified, 'nonaffine oracle regression/shift/reflection')
    print('ORACLE VERIFIED', summary['max_model_error_N_m3'], summary['BULK_MODEL_DIFFERENCE'], flush=True)


def decomposition_stage(cfg):
    require(load('BULK_ORACLE_VALIDATION.json')['VERIFIED'], 'oracle required first')
    field = QuadraticAxial(); h = cfg['h_m']; Fref = field.scales(h)['F_ref_N']
    # Preserve the previous round's exact floating construction L=(N-1)*h.
    g, f = block_family(cfg['N'], (cfg['N']-1)*h, cfg['delta_m']); mat = intact_material()
    # Fresh affine finite-domain evaluation, compared to the frozen previous B1.
    aff = QuadraticAxial(beta=0.); st = evaluate(f, aff.position(f.X), mat)
    raw = load('BULK_ORACLE_RAW.json'); qa = np.array(next(r['q_bulk_N_m3'] for r in raw if r['label']=='affine' and r['level']==9))
    with np.load(OLD/'NORMAL_AFFINE.npz') as old:
        require(np.array_equal(old['X'], f.X), 'old/new physical grid mismatch')
        affine_load = qa*f.volume[:, None]-st.force
        err = float(np.linalg.norm(affine_load-old['load'])/np.linalg.norm(old['load']))
        old_error = normmax(-st.force-old['load'])/Fref
    save('AFFINE_LIMIT_RECOVERY.json', dict(relative_L2_error=err, old_force_difference_scaled=old_error, q_bulk_N_m3=qa,
                                          gate=cfg['gate'], PASS=err<=cfg['gate']))
    require(err<=cfg['gate'], 'affine B1 recovery')
    del st; gc.collect()
    target = field.position(f.X); st = evaluate(f, target, mat)
    zs, iz = np.unique(f.X[:,2], return_inverse=True)
    q = np.array([r['q_bulk_N_m3'] for r in raw if r['label']=='primary'])[iz]
    b = field.response(f.X[:,2])['b']; arrays = localize(st.force, q, b, f.volume)
    closure = normmax(arrays['f_FD']+arrays['f_BC_localized']+arrays['f_body_local']-arrays['f_model'])/Fref
    naive = normmax(arrays['f_BC_naive']-arrays['f_BC_localized']+arrays['f_bulk_target'])/Fref
    nl = normmax(st.force+arrays['f_external_NL'])/Fref
    analyticF = field.response(f.X[:,2])['F']
    target_metrics = dict(energy_J=st.energy, Fbar_vs_analytic_max_scaled=float(np.linalg.norm(st.Fbar-analyticF,axis=(1,2)).max()/field.scales(h)['gradient_scale']))
    # Check the nonaffine selective adaptation against the full evaluator, without a new method.
    selective = ManufacturedStencil(field, 19); selective_checks=[]
    for offset in [[0,0,0],[-5,0,0],[0,0,8]]:
        pid = selective.id_at_offset(offset); fs = selective.force(pid)
        selective_checks.append(dict(offset=offset, ID=pid, force_N=fs,
                                     difference_scaled=float(np.linalg.norm(fs-st.force[pid])/Fref)))
    save('SELECTIVE_VS_FULL_NONAFFINE.json', selective_checks)
    require(max(r['difference_scaled'] for r in selective_checks)<=cfg['gate'], 'selective vs full nonaffine')
    np.savez_compressed(OUT/'NODEWISE_RESIDUAL_DECOMPOSITION.npz', X=f.X, volume=f.volume,
                        x_target=target, q_bulk=q, b_local=b, e_model=q+b, Fbar_target=st.Fbar,
                        F_local=analyticF, distance_faces=g.face_distances(f.X), d_min_over_delta=g.distance(f.X)/f.horizon,
                        d_min_over_h=g.distance(f.X)/h, **arrays)
    save('NODEWISE_IDENTITIES.json', dict(closure_scaled=closure, naive_difference_scaled=naive, NL_target_residual_scaled=nl,
                                         target_metrics=target_metrics, PASS=max(closure,naive,nl)<=cfg['gate']))
    require(max(closure,naive,nl)<=cfg['gate'], 'nodewise residual decomposition')
    force,moment,work = balances(f.X,arrays,field.continuum(),Fref,field.length)
    save('FORCE_BALANCE.json', force); save('MOMENT_BALANCE.json', moment); csvsave('VIRTUAL_WORK_DECOMPOSITION.csv', work)
    rows=[]; quantiles={}
    for name in ['localized','naive']:
        profile, qs = support_profile(g,f,arrays['f_BC_'+name],h)
        for qv in qs.values(): qv['depth_um']=qv['depth_m']*1e6
        quantiles[name]=qs
        rows.extend([dict(correction=name,**r) for r in profile])
    csvsave('BOUNDARY_LOCALIZATION.csv',rows); save('SUPPORT_QUANTILES.json',quantiles)
    csvsave('NODEWISE_RESIDUAL_SUMMARY.csv',[dict(level=k,Z_m=z,nodes=int((iz==k).sum()),
                **{name+'_z_mean_N':float(a[iz==k,2].mean()) for name,a in arrays.items()},
                **{name+'_norm_max_N':normmax(a[iz==k]) for name,a in arrays.items()}) for k,z in enumerate(zs)])
    print('DECOMPOSITION VERIFIED',closure,'AFFINE',err,'TARGET_RECONSTRUCTION',target_metrics,flush=True)
    del st,selective,f;gc.collect()


def n27_stage(cfg):
    require(load('NODEWISE_IDENTITIES.json')['PASS'], 'nodewise required first')
    field=QuadraticAxial();s=ManufacturedStencil(field,27);Fref=field.scales(cfg['h_m'])['F_ref_N']
    # Independent enlarged/shifted fully closed scaffold for the center comparison.
    reference = oracle_job(dict(label='independent_N29',Z=0.,n=29,offset=[1,-1,1]))
    q=np.asarray(reference['q_bulk_N_m3']); rows=[]; closures=[]
    for offset in cfg['N27_points_offsets']:
        pid=s.id_at_offset(offset);f=s.force(pid);bc=q*s.V[pid]-f;cl=s.closure(pid);closures.append(cl)
        depth=int(np.minimum(s.ijk[pid],26-s.ijk[pid]).min())
        rows.append(dict(offset=offset,point_X_m=s.X[pid],depth_h=depth,depth_over_delta=depth*s.h/s.delta,
                         force_FD_N=f,bulk_target_N=q*s.V[pid],localized_N=bc,localized_norm_N=float(np.linalg.norm(bc)),
                         localized_scaled=float(np.linalg.norm(bc)/Fref),naive_norm_N=float(np.linalg.norm(f)),
                         naive_scaled=float(np.linalg.norm(f)/Fref),fully_closed=cl['FULL_FORCE_STENCIL_CLOSED']))
        print('N27',offset,rows[-1]['localized_scaled'],cl['FULL_FORCE_STENCIL_CLOSED'],flush=True)
    save('N27_DEPTH_CHECK.json',dict(points=rows,independent_reference=reference,closure=closures,
                                   PASS=rows[-1]['fully_closed'] and rows[-1]['localized_scaled']<=cfg['gate']))
    require(rows[-1]['fully_closed'] and rows[-1]['localized_scaled']<=cfg['gate'], 'N27 closed center localized correction')


def solve_stage(cfg,kind):
    require(load('N27_DEPTH_CHECK.json')['PASS'], 'N27 gate required before solves')
    if kind=='classical': require(load('NL_MANUFACTURED_RECOVERY.json')['PASS'],'nonlocal recovery required before classical solve')
    field=QuadraticAxial();h=cfg['h_m'];sc=field.scales(h);g,f=block_family(19,18*h,cfg['delta_m'])
    with np.load(OUT/'NODEWISE_RESIDUAL_DECOMPOSITION.npz') as saved:
        target=saved['x_target'].copy();Ftarget=saved['Fbar_target'].copy();analytic=saved['F_local'].copy()
        external=saved['f_external_NL' if kind=='NL' else 'f_external_classical'].copy()
    fixed=g.dirichlet_mask(f.X);q=(f.X-g.lower)/g.lengths
    perturb=.02*cfg['epsilon']*cfg['L_m']*np.prod(np.sin(np.pi*q),axis=1)
    initial=target+perturb[:,None]*np.array([1.,-.4,.7]);initial[fixed]=target[fixed]
    np.savez_compressed(OUT/f'{kind}_SOLVE_INPUT.npz',initial=initial,external=external,fixed=fixed)
    print('START SOLVE',kind,'free DOFs',3*int((~fixed).sum()),flush=True)
    x,info=solve(f,intact_material(),target,fixed,initial,external,h,sc['F_ref_N'],cfg)
    st=evaluate(f,x,intact_material());res=st.force+external
    solver_scale=max(float(np.linalg.norm(external[~fixed])),sc['F_ref_N'])
    fr=float(np.linalg.norm(res[~fixed])/solver_scale)
    ferr=np.linalg.norm(st.Fbar-Ftarget,axis=(1,2))/sc['gradient_scale']
    fanalytic=np.linalg.norm(st.Fbar-analytic,axis=(1,2))/sc['gradient_scale']
    uerr=np.linalg.norm(x-target,axis=1)
    target_energy=load('NODEWISE_IDENTITIES.json')['target_metrics']['energy_J']
    result=dict(kind=kind,solve=info,force_residual=fr,solver_force_scale_N=solver_scale,
                force_max_free_local_scaled=normmax(res[~fixed])/sc['F_ref_N'],
                F_recovery_error_max=float(ferr.max()),F_analytic_error_max=float(fanalytic.max()),
                imposed_gradient_scale=sc['gradient_scale'],displacement_error_max_m=float(uerr.max()),
                displacement_relative_L2=float(np.linalg.norm(x-target)/np.linalg.norm(target-f.X)),
                initial_perturbation_max_m=normmax(initial-target),energy_J=float(st.energy),target_energy_J=target_energy,
                energy_deviation_J=float(st.energy-target_energy),energy_relative_deviation=float((st.energy-target_energy)/target_energy),
                fixed_nodes=int(fixed.sum()),fixed_displacement_error_m=float(abs(x[fixed]-target[fixed]).max()),
                PASS=bool(info['converged'] and fr<=cfg['mechanical_residual_gate'] and (kind!='NL' or ferr.max()<=cfg['gradient_recovery_gate'])))
    save('NL_MANUFACTURED_RECOVERY.json' if kind=='NL' else 'CLASSICAL_MANUFACTURED_RESULT.json',result)
    np.savez_compressed(OUT/f'{kind}_RECOVERED_STATE.npz',X=f.X,x=x,volume=f.volume,initial=initial,target=target,
                        external=external,force_internal=st.force,reaction=-res*fixed[:,None],fixed=fixed,Fbar=st.Fbar,
                        F_recovery_error=ferr,F_analytic_error=fanalytic,u_error=uerr)
    depth=np.rint(g.distance(f.X)/h).astype(int)
    csvsave(f'{kind}_RECOVERY_PROFILE.csv',[dict(depth_h=int(k),depth_over_delta=k*h/f.horizon,nodes=int((depth==k).sum()),
                F_recovery_max=float(ferr[depth==k].max()),F_recovery_RMS=float(np.sqrt(np.mean(ferr[depth==k]**2))),
                F_analytic_max=float(fanalytic[depth==k].max()),u_error_max_m=float(uerr[depth==k].max())) for k in np.unique(depth)])
    print('DONE SOLVE',kind,'force',fr,'max F recovery',ferr.max(),'max F analytic',fanalytic.max(),flush=True)
    if kind=='NL':require(result['PASS'],'nonlocal manufactured perturbed-start recovery')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['oracle','decomposition','n27','NL','classical'],required=True)
    args=ap.parse_args();cfg=load('PREDECLARED_CONFIG.json');require(integrity()['PASS'],'protected source/input integrity')
    require(psutil.virtual_memory().available>=.2*psutil.virtual_memory().total,'RAM headroom before stage')
    started=time.monotonic();stop=threading.Event();p=psutil.Process()
    stats=dict(stage=args.stage,CPU_workers=cfg['oracle_workers'] if args.stage=='oracle' else 1,BLAS_threads=1,
               GPU='GPU_NOT_USED_BY_DESIGN',peak_process_tree_RSS_bytes=0,minimum_available_RAM_bytes=psutil.virtual_memory().available)
    def monitor():
        while not stop.is_set():
            rss=0
            for proc in [p]+p.children(recursive=True):
                try:rss+=proc.memory_info().rss
                except psutil.Error:pass
            stats['peak_process_tree_RSS_bytes']=max(stats['peak_process_tree_RSS_bytes'],rss)
            stats['minimum_available_RAM_bytes']=min(stats['minimum_available_RAM_bytes'],psutil.virtual_memory().available)
            stop.wait(.1)
    t=threading.Thread(target=monitor,daemon=True);t.start()
    try:
        if args.stage=='oracle':oracle_stage(cfg)
        elif args.stage=='decomposition':decomposition_stage(cfg)
        elif args.stage=='n27':n27_stage(cfg)
        else:solve_stage(cfg,args.stage)
    finally:
        stop.set();t.join();stats['wall_seconds']=time.monotonic()-started
        stats['peak_main_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        stats['minimum_available_RAM_fraction']=stats['minimum_available_RAM_bytes']/psutil.virtual_memory().total
        save(f'RESOURCES_{args.stage}.json',stats)
        save('HISTORICAL_INTEGRITY.json',integrity())
    require(stats['minimum_available_RAM_fraction']>=.2,'RAM headroom during stage')


if __name__=='__main__':main()
