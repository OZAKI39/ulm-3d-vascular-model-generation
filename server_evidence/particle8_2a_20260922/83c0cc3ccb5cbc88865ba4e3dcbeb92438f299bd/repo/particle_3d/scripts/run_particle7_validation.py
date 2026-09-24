#!/usr/bin/env python3
from pathlib import Path
import sys,argparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle7_validation import *
p=argparse.ArgumentParser(); p.add_argument('--real-only',action='store_true'); p.add_argument('--isolated-mb-only',action='store_true'); p.add_argument('--output',default=str(REPO/'particle_3d/reports/particle7/data'));  a=p.parse_args()
out=Path(a.output).resolve(); out.mkdir(parents=True,exist_ok=True)
if a.isolated_mb_only:
 r=isolated_real_mb_transport(out); print('ISOLATED_MB',r['status'],r['displacement_m'],flush=True)
elif a.real_only:
 r=real_smoke(out); print('REAL',r['accounting'],r['rejection_counts'],flush=True)
else:
 audit,samplers=inlet_audit_validation(out); scheduler_validation(out,audit['INLET']['positive_Q_m3_s']); print('SCHEDULERS',flush=True)
 sampler_validation(out,samplers['INLET']); print('SAMPLER_ORIENTATION',flush=True)
 admission_validation(out); print('ADMISSION',flush=True)
 lifecycle_restart_validation(out); print('LIFECYCLE_RESTART',flush=True)
