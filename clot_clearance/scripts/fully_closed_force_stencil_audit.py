"""Focused fully closed BAQP force-stencil audit; no historical/full-suite QA."""
import csv
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import threading
import time
import numpy as np
from pd_clot.pfpd_candidate.force_stencil_audit import (
    H,DELTA,EPSILON,FORCE_TOL,COEFFICIENT_TOL,RECONSTRUCTION_TOL,
    configuration,bulk_oracle,AffineStencil,validate_selective)

ROOT=Path(__file__).resolve().parents[1]
NAME='pfpd_baqp_fully_closed_force_stencil_audit_v1'
OUT=ROOT/'verification'/NAME
RESULTS=ROOT/'results'/NAME


def encode(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value))


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,default=encode,indent=2,allow_nan=False)+'\n')


def table(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


class Resources:
    def __init__(self):
        self.start=time.monotonic();self.stop=threading.Event();self.peak=0;self.min_available=10**30
        self.total=self.memory()['MemTotal'];self.thread=threading.Thread(target=self.monitor,daemon=True)
        if self.memory()['MemAvailable']<.2*self.total:raise RuntimeError('Insufficient 20% RAM headroom')
        self.thread.start()

    def memory(self):
        return {s.split(':')[0]:int(s.split()[1])*1024 for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith(('MemTotal:','MemAvailable:'))}

    def tree_rss(self,pid):
        try:
            status=Path(f'/proc/{pid}/status').read_text().splitlines()
            rss=next((int(s.split()[1])*1024 for s in status if s.startswith('VmRSS:')),0)
            children=Path(f'/proc/{pid}/task/{pid}/children').read_text().split()
            return rss+sum(self.tree_rss(int(p)) for p in children)
        except (FileNotFoundError,ProcessLookupError):return 0

    def monitor(self):
        while not self.stop.is_set():
            self.peak=max(self.peak,self.tree_rss(os.getpid()))
            self.min_available=min(self.min_available,self.memory()['MemAvailable'])
            self.stop.wait(.05)

    def finish(self):
        self.stop.set();self.thread.join()
        return dict(CPU_workers=1,BLAS_threads=1,GPU='NOT_USED_BY_DESIGN',
                    wall_seconds=time.monotonic()-self.start,process_tree_peak_RSS_sampled_bytes=self.peak,
                    process_peak_RSS_rusage_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                    total_RAM_bytes=self.total,minimum_available_RAM_sampled_bytes=self.min_available,
                    minimum_available_fraction=self.min_available/self.total,
                    RAM_headroom_20percent_maintained=self.min_available>=.2*self.total,
                    sample_interval_seconds=.05,full_N19_N27_bond_states_allocated=False)


def observation_offsets(n):
    if n!=27:return [('center',(0,0,0)),('near_center_bottom',(0,0,-1)),('near_center_top',(0,0,1))]
    m=n//2;depths=[('near_1delta',int(round(DELTA/H))),('closed_1hop',int(np.ceil(DELTA/H))),
                  ('near_2delta',int(round(2*DELTA/H))),('closed_2hop',int(np.ceil(2*DELTA/H))),
                  ('near_3delta',int(round(3*DELTA/H)))]
    out=[('center',(0,0,0))]
    for label,depth in depths:
        out.extend([(label+'_bottom',(0,0,depth-m)),(label+'_top',(0,0,m-depth))])
    return out


def case(n):
    start=time.monotonic();grid=AffineStencil(n);observations=[]
    for label,offset in observation_offsets(n):
        pid=grid.id_at_offset(offset)
        result=grid.evaluate_point(pid);result['label']=label
        result['closure']=grid.closure(pid)
        result['coefficients']=grid.coefficients(pid)
        observations.append(result)
    center=observations[0]
    symmetry=[]
    for bottom in [r for r in observations if r['label'].endswith('_bottom')]:
        top=next(r for r in observations if r['label']==bottom['label'].replace('_bottom','_top'))
        defect=top['force_N']+bottom['force_N']
        symmetry.append(dict(bottom=bottom['label'],top=top['label'],pair_defect_N=defect,
                             pair_defect_norm_N=float(np.linalg.norm(defect)),
                             pair_defect_local_normalized=float(np.linalg.norm(defect)/grid.scale)))
    diagnostic=grid.stress_roundoff_diagnostic(center['ID'])
    return dict(n=n,L_m=grid.L,particles=n**3,sum_volume_m3=float(grid.V.sum()),
                volume_relative_error=float(abs(grid.V.sum()/grid.L**3-1)),observations=observations,
                symmetry_pairs=symmetry,center_stress_roundoff_diagnostic=diagnostic,
                selective_cache_counts=grid.cache_counts(),wall_seconds=time.monotonic()-start)


def status_for(cases,oracle,validation):
    center={c['n']:c['observations'][0] for c in cases};n27=center[27]
    errors=n27['coefficients']['max_relative_error']
    closure=n27['closure']['FULL_FORCE_STENCIL_CLOSED']
    force_ok=n27['force_local_normalized']<=FORCE_TOL
    coefficients_ok=max(errors.values())<=COEFFICIENT_TOL
    if closure and force_ok and coefficients_ok:
        patch='PASS';root='FINITE_DOMAIN_MULTI_HOP_BOUNDARY_INFLUENCE'
    elif closure and not force_ok:
        patch='FAIL';root='SELECTIVE_EVALUATOR_INCONSISTENCY' if oracle['PASS'] else 'BAQP_BULK_AFFINE_FORCE_PATCH_FAILURE'
    elif not closure:
        patch='NOT_ESTABLISHED';root='FORCE_STENCIL_DEPTH_UNDERESTIMATED'
    else:
        patch='NOT_ESTABLISHED';root='INCONCLUSIVE'
    off=[r for c in cases for r in c['observations'] if r['label']!='center']
    has_z=any(np.linalg.norm(r['channels']['original']['Z_correction'])/oracle['local_force_scale_N']>FORCE_TOL for r in off)
    has_excess=any(np.linalg.norm(r['channels']['normalization_excess']['force'])/oracle['local_force_scale_N']>FORCE_TOL for r in off)
    secondary='BOTH' if has_z and has_excess else 'Z_MEDIATED_PROPAGATION' if has_z else 'OMEGA0_NORMALIZATION_PROPAGATION' if has_excess else 'NONE'
    status=dict(INFINITE_LATTICE_BULK_ORACLE='PASS' if oracle['PASS'] else 'FAIL',
                SELECTIVE_EVALUATOR='VERIFIED' if all(v['PASS'] for v in validation) else 'NOT_VERIFIED')
    for n,c in center.items():
        status[f'N{n}_CLOSURE']='/'.join(f"D{r['depth']}:{'CLOSED' if r['closed_with_bulk_volumes'] else 'OPEN'}" for r in c['closure']['levels'][1:])
    status.update(FULL_FORCE_STENCIL_CLOSED='YES' if closure else 'NO',ACTUAL_REQUIRED_GRAPH_DEPTH=3,
                  MINIMUM_OBSERVED_CLOSURE_DEPTH_FOR_BULK_FORCE='NOT_ESTABLISHED',
                  center_force_only_first_observed_closed_hop=1,
                  minimum_depth_reason='Every centered cube passes the force-only gate by inversion symmetry. Off-center two-hop points do not all pass. A minimum generally sufficient closure cannot be inferred from center cancellation or a few passing points.',
                  N27_CENTER_FORCE_N=n27['force_norm_N'],N27_CENTER_FORCE_VECTOR_N=n27['force_N'],
                  N27_CENTER_FORCE_LOCAL_NORMALIZED=n27['force_local_normalized'],
                  N27_CENTER_FORCE_ONLY_GATE='PASS' if force_ok else 'FAIL',
                  N27_COEFFICIENT_GATE='PASS' if coefficients_ok else 'FAIL',
                  BULK_AFFINE_FORCE_PATCH_TEST=patch,OMEGA0_BULK_MATCH=errors['omega0'],K_BULK_MATCH=errors['K'],Z_BULK_MATCH=errors['Z'],
                  ROOT_CAUSE_PRIMARY=root,ROOT_CAUSE_SECONDARY=secondary,
                  CURRENT_C5_C7_ACTIVATION_RESULTS='UNCHANGED',BOUNDARY_FIX='NOT_IMPLEMENTED',FRACTURE='NOT_RUN',C7_STABILITY='NOT_RUN')
    return status


def export_tables(cases,scale):
    points=[];channels=[];centers=[];coeff=[]
    for c in cases:
        for r in c['observations']:
            f=r['force_N'];lev=r['closure']['levels'];co=r['coefficients']['max_relative_error']
            row=dict(domain=f"N{c['n']}",point=r['label'],ID=r['ID'],distance_over_delta=r['boundary_distance_over_delta'],
                     distance_h=r['boundary_distance_cells'],fx_N=f[0],fy_N=f[1],fz_N=f[2],norm_N=r['force_norm_N'],
                     norm_local_scaled=r['force_local_normalized'],abs_fz_local_scaled=r['abs_fz_local_normalized'],
                     D1_closed=lev[1]['closed_with_bulk_volumes'],D2_closed=lev[2]['closed_with_bulk_volumes'],D3_closed=lev[3]['closed_with_bulk_volumes'],
                     full_stencil_closed=r['closure']['FULL_FORCE_STENCIL_CLOSED'])
            points.append(row)
            if r['label']=='center':centers.append(row)
            coeff.append(dict(domain=row['domain'],point=r['label'],**co))
            for weight,data in r['channels'].items():
                for component,vec in data.items():
                    channels.append(dict(domain=row['domain'],point=r['label'],weight=weight,component=component,
                                         fx_N=vec[0],fy_N=vec[1],fz_N=vec[2],norm_N=float(np.linalg.norm(vec)),
                                         norm_local_scaled=float(np.linalg.norm(vec)/scale),fz_local_scaled=float(vec[2]/scale)))
    for name,rows in [('CENTER_AFFINE_FORCE.csv',centers),('OBSERVATION_POINTS.csv',points),('FORCE_DECOMPOSITION.csv',channels),('COEFFICIENT_CONVERGENCE.csv',coeff)]:
        table(RESULTS/'tables'/name,rows)


def terminal_summary(summary):
    c=summary['configuration'];o=summary['oracle'];s=summary['status'];r=summary['resources'];git=summary['git']
    lines=[f"branch / HEAD: {git['branch'].strip()} / {git['HEAD'].strip()}",f'h: {H:.12g} m',f'delta: {DELTA:.12g} m',f'delta/h: {DELTA/H:g}',f'epsilon: {EPSILON:g}',
           f"Pzz: {o['Pzz']:.16g} Pa",f"local force scale: {o['local_force_scale_N']:.16e} N",
           f"infinite lattice: offsets={o['neighbor_offset_count']}; force={o['center_force_N']}; normalized={o['center_force_local_normalized']:.8e}"]
    for case_data in summary['cases']:
        point=case_data['observations'][0];n=case_data['n']
        lines.append(f"N{n}: d/delta={point['boundary_distance_over_delta']:.9g}; {s[f'N{n}_CLOSURE']}; center force/local={point['force_local_normalized']:.8e}")
    lines.extend([f"coefficient bulk errors: omega0={s['OMEGA0_BULK_MATCH']:.8e}; K={s['K_BULK_MATCH']:.8e}; Z={s['Z_BULK_MATCH']:.8e}",
                  f"minimum closure depth: {s['MINIMUM_OBSERVED_CLOSURE_DEPTH_FOR_BULK_FORCE']}",f"BULK_AFFINE_FORCE_PATCH_TEST: {s['BULK_AFFINE_FORCE_PATCH_TEST']}",
                  f"ROOT_CAUSE_PRIMARY: {s['ROOT_CAUSE_PRIMARY']}",f"ROOT_CAUSE_SECONDARY: {s['ROOT_CAUSE_SECONDARY']}",
                  f"CPU workers: {r['CPU_workers']}",f"wall time: {r['wall_seconds']:.3f} s (recorded audit computation, targeted tests, figure exports)",
                  f"peak RAM: {r['process_tree_peak_RSS_sampled_bytes']/2**30:.3f} GiB (sampled process tree)",f"GPU: {r['GPU']}",
                  f"REVIEW: {OUT/'REVIEW_FULLY_CLOSED_FORCE_STENCIL_V1.md'}",f"figures: {RESULTS/'figures'}"])
    return '\n'.join(lines)+'\n'


def main():
    resources=Resources();OUT.mkdir(parents=True,exist_ok=True)
    config=configuration();existing=OUT/'PREDECLARED_CONFIG.json'
    if existing.exists():
        assert json.loads(existing.read_text())==config,'Do not alter predeclared controls'
    else:save(existing,config)
    oracle=bulk_oracle();save(OUT/'BULK_LATTICE_ORACLE.json',oracle)
    if not oracle['PASS']:
        save(OUT/'STOP_GATE.json',dict(status='INFINITE_LATTICE_BULK_FORCE_NONZERO'));raise SystemExit(1)
    validation=[]
    for n in (11,13):
        v=validate_selective(n);validation.append(v);save(OUT/f'SELECTIVE_VALIDATION_N{n}.json',v)
        if not v['PASS']:
            save(OUT/'STOP_GATE.json',dict(status='SELECTIVE_EVALUATOR_INCONSISTENCY',n=n));raise SystemExit(1)
    data=[case(n) for n in (11,19,27)]
    for c in data:save(OUT/f"N{c['n']}.json",c)
    save(OUT/'DEPENDENCY_CLOSURE.json',dict(cases=[dict(n=c['n'],points=[dict(label=p['label'],**p['closure']) for p in c['observations']]) for c in data]))
    save(OUT/'COEFFICIENT_CONVERGENCE.json',dict(cases=[dict(n=c['n'],points=[dict(label=p['label'],**p['coefficients']) for p in c['observations']]) for c in data]))
    zero=AffineStencil(11,epsilon=0.);zero_rows=[zero.evaluate_point(zero.id_at_offset(v)) for v in ((0,0,0),(0,0,1),(1,-1,-1))]
    noise=dict(local_force_scale_N=oracle['local_force_scale_N'],zero_strain_normalization='Fixed nonzero loaded-state scale, never divide by Pzz(0)=0',
               bulk_oracle_force_N=oracle['center_force_N'],bulk_oracle_normalized=oracle['center_force_local_normalized'],
               zero_strain_max_force_N=max(r['force_norm_N'] for r in zero_rows),zero_strain_max_local_normalized=max(r['force_local_normalized'] for r in zero_rows),
               symmetry_pairs=[dict(n=c['n'],pairs=c['symmetry_pairs']) for c in data],
               N27_stress_reconstruction=data[-1]['center_stress_roundoff_diagnostic'])
    save(OUT/'NOISE_BASELINE.json',noise)
    status=status_for(data,oracle,validation);save(OUT/'STATUS.json',status)
    export_tables(data,oracle['local_force_scale_N'])
    # Only the seven requested tests. No repository-wide or historical suite.
    t=time.monotonic()
    run=subprocess.run([sys.executable,'-B','-m','pytest','tests/test_fully_closed_force_stencil.py','-q','--junitxml='+str(OUT/'TARGETED_TESTS.xml')],
                       cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (OUT/'TARGETED_TESTS.log').write_text(run.stdout)
    test_info=dict(command='pytest tests/test_fully_closed_force_stencil.py -q',returncode=run.returncode,wall_seconds=time.monotonic()-t)
    save(OUT/'TARGETED_TESTS.json',test_info)
    if run.returncode:raise RuntimeError('Targeted tests failed; see TARGETED_TESTS.log')
    summary=dict(configuration=config,git=json.loads((OUT/'STARTING_GIT_STATE.json').read_text()),oracle=oracle,
                 selective_validation=validation,cases=data,noise=noise,status=status,targeted_tests=test_info)
    save(OUT/'SCIENTIFIC_RESULTS.json',summary)
    save(OUT/'COMPUTE_RESOURCES.json',dict(CPU_workers=1,wall_seconds=time.monotonic()-resources.start,
         process_tree_peak_RSS_sampled_bytes=resources.peak,minimum_available_RAM_sampled_bytes=resources.min_available,
         total_RAM_bytes=resources.total))
    from scripts.plot_closed_force_stencil import make_figures,write_review
    summary['figures']=make_figures(summary,RESULTS/'figures')
    summary['resources']=resources.finish()
    if not summary['resources']['RAM_headroom_20percent_maintained']:raise RuntimeError('RAM headroom requirement not met')
    save(OUT/'RESOURCES.json',summary['resources']);save(OUT/'SUMMARY.json',summary)
    write_review(summary,OUT/'REVIEW_FULLY_CLOSED_FORCE_STENCIL_V1.md',RESULTS/'figures')
    text=terminal_summary(summary);(OUT/'TERMINAL_SUMMARY.txt').write_text(text);print(text,end='')


if __name__=='__main__':main()
