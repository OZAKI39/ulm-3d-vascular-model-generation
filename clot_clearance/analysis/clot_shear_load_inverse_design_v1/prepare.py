from common import *
from shear_field_library import freeze_library
from current_surface import freeze_surface
from multifragment_shear_objective import freeze
import shutil,subprocess,sys,importlib.metadata,datetime

assert not (O/'BASELINE_MANIFEST.json').exists(),'Do not overwrite a prior preparation'
assert sha(PREV/'DELIVERY_MANIFEST.json')=='ba855a0de4351e14d32d96d31a3674f5c8deccc85eafc1195eb751e6f2414317'
assert sha(LIB)=='465f8144c6da5c9de3aec4137b9c325f53c26737b199c536fb4a20845511ef77'
(O/'provenance').mkdir(parents=True,exist_ok=True)
for name in ['reference.npz','reference.native','body.data']:
    shutil.copy2(PREV/'provenance'/name,O/'provenance'/name)
dependencies=[LIB,PREV/'DELIVERY_MANIFEST.json',PREV/'SHEAR_TRACTION_FIELD.npz',
              R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',
              R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',PA/'native.py',
              *(O/'provenance').iterdir()]
baseline=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files={str(p):dict(sha256=sha(p),bytes=p.stat().st_size) for p in dependencies})
dump(O/'BASELINE_MANIFEST.json',baseline)
ref=load_npz(O/'provenance/reference.npz');freeze_library();freeze_surface(ref);freeze()
dump(O/'MATERIAL_AND_SPECIMEN_FREEZE.json',dict(material_nodes=int(ref['clot_count']),anchor_nodes=len(ref['X'])-int(ref['clot_count']),original_bonds=len(ref['bonds']),initial_precuts=int(ref['initial_cut'].sum()),G_Pa=G,K_Pa=K,rho_kg_m3=RHO,dx_m=DX,horizon_m=HORIZON,bulk_s0=S0[1],bulk_sc=SC[1],interface_s0=S0[2],interface_sc=SC[2],anchor_s0=S0[3],anchor_sc=SC[3],heterogeneity='NONE_ORIGINAL_UNIFORM',law='Unchanged native variational LPS, kappa/d irreversible, old precrack preserved',contact='Unchanged bpm/peri, contact_distance_factor=1.0, nonbonded contact does not define connectivity',status='EXPLORATORY_UNCALIBRATED'))
versions={p:importlib.metadata.version(p) for p in ['numpy','scipy','matplotlib','mpi4py','networkx']}
dump(O/'ENVIRONMENT.json',dict(python=sys.version,executable=sys.executable,packages=versions,rank_count=1,threads=1,new_installs=[]))
print('PREPARED',TASK,RUN,len(ref['X']),len(ref['bonds']))
