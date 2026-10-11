"""Prepare a NEW explicitly named replay; never overwrite the accepted run.
Run using the existing activated server environment. Subsequent solvers require
individual task-specific Supervisor jobs and the same conditional entry gates.
"""
import os,argparse
parser=argparse.ArgumentParser();parser.add_argument('--run-id',required=True);args=parser.parse_args()
assert args.run_id!='20261011T154854Z' and '/' not in args.run_id and '..' not in args.run_id
os.environ['Y_CLOT_DUAL_RUN_ID']=args.run_id
from analysis.dual_track_v1.common import *
for p in [A,B,COMMON]:p.mkdir(parents=True,exist_ok=False)
from analysis.dual_track_v1.geometry_b import surface
surface()
# Native-compatible tag representation produced directly by the corrected source.
import shutil,pyvista as pv
source=B/'clot_present_lumen.vtp';target=B/'clot_present_lumen_int32.vtp';s=pv.read(source);assert s['ModelFaceID'].dtype==np.int32;shutil.copy2(source,target)
save(B/'GEOMETRY_REPRESENTATION_ADDENDUM.json',dict(old_sha256=sha(source),new_sha256=sha(target),values_identical=True,formal_CFD_attempts_at_change=0,reason='Corrected replay directly emits the native int32 tag format. No shape edits.'))
from analysis.dual_track_v1.design import main
main()
print('Prepared only. Run all subsequent commands with Y_CLOT_DUAL_RUN_ID='+args.run_id)
