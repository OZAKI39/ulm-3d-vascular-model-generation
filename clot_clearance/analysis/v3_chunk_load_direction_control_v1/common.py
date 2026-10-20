"""Isolated paths; unchanged validated V3 constants and numerical implementations."""
from pathlib import Path
import importlib.util,sys,os
R=Path(os.environ.get('Y_CLOT_PROJECT_ROOT','/workspace/y_clot/project'))
HIST_A=R/'analysis/v3_softening_pilot_v1'
def module_from_file(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
_old=module_from_file('direction_frozen_common',HIST_A/'common.py')
for _k,_v in vars(_old).items():
 if not _k.startswith('_'):globals()[_k]=_v
R=Path(os.environ.get('Y_CLOT_PROJECT_ROOT','/workspace/y_clot/project'));HIST_A=R/'analysis/v3_softening_pilot_v1'
TASK='v3_chunk_fragmentation_load_direction_control_v1';RUN='20261016T105642Z'
A=Path(__file__).resolve().parent;O=R/'verification'/TASK/RUN;F=R/'results'/TASK/RUN
B=R/'verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z'
PREV=R/'verification/v3_chunk_fragmentation_mechanism_v1/20261015T174630Z';PREV_A=R/'analysis/v3_chunk_fragmentation_mechanism_v1'
LIB=B/'build/liblammps.so';EARLY=.0008
SURFACE=R/'results/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/paraview'
sys.path.insert(1,str(HIST_A));sys.path.insert(2,str(PREV_A))
def load_npz(p):
 with np.load(p) as z:return {k:z[k] for k in z.files}
