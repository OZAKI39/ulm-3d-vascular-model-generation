#!/usr/bin/env python3
"""Stage B then C: direct point basin for every retained P8.1 proposal."""
from pathlib import Path
from types import SimpleNamespace
import argparse,json,multiprocessing,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_diagnostics import prepare_p81_proposals,summarize_admission,stop_audit
from particle_3d.particle82_tracers import run,init_environment
from particle_3d.particle82_provenance import atomic_json,require_remote
from particle_3d.injection_admission import FiniteSizeAdmission

ENV=None

def recheck(meta):
    checker=FiniteSizeAdmission(ENV.sampler,SimpleNamespace(),wall=ENV.wall,field=ENV.field)
    errors=[]
    for row in meta['admission_candidates']:
        _,status,_=checker.check(meta['birth_metadata'],row['position_m'],{})
        if status!=row['status']:errors.append(dict(draw=row['draw'],original=row['status'],rechecked=status))
    return dict(particle_id=meta['particle_id'],checked=len(meta['admission_candidates']),mismatches=errors)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--previous',required=True);p.add_argument('--output-dir',required=True)
    p.add_argument('--host-provenance',required=True);p.add_argument('--point-basin-summary',required=True)
    p.add_argument('--workers',type=int,default=8);a=p.parse_args()
    host=json.loads(Path(a.host_provenance).read_text());require_remote(host,host['hostname'])
    point=json.loads(Path(a.point_basin_summary).read_text())
    if point['count']<100000 or not point['three_outlets_observed']:raise ValueError('Blocking point-tracer three-outlet gate')
    out=Path(a.output_dir);positions=prepare_p81_proposals(a.previous,out)
    ENV=init_environment()
    meta=[json.loads(f.read_text()) for f in sorted((Path(a.previous)/'trajectories').glob('*.json'))]
    # Re-execute the original check at the exact saved proposals/radii. Do not
    # move points or rerun the RNG to find more acceptable positions.
    with multiprocessing.get_context('fork').Pool(a.workers) as pool:checks=list(pool.imap_unordered(recheck,meta,chunksize=1))
    checks.sort(key=lambda r:r['particle_id'])
    audit=dict(all_exact=all(not r['mismatches'] for r in checks),checked_proposals=sum(r['checked'] for r in checks),rows=checks)
    atomic_json(out/'ORIGINAL_ADMISSION_RECHECK.json',audit)
    if not audit['all_exact']:raise ValueError('Original admission recheck mismatch')
    run(positions,out/'point_paths',host,workers=a.workers,shard_size=256,config=dict(step_m=1e-7,error=1e-13,horizon_m=2e-3))
    admission=summarize_admission(a.previous,out,out/'point_paths')
    stop=stop_audit(a.previous,admission,out/'stop_audit')
    print(json.dumps(dict(admission=admission['by_basin'],stops=stop['by_basin'])),flush=True)
