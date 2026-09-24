#!/usr/bin/env python3
"""Fresh render windows must replay the same saved active state identically."""
from pathlib import Path
import argparse,hashlib,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from particle_3d.particle82_results import Scene
from particle_3d.particle82_visuals import compose,schedule
from particle_3d.particle81_visuals import VesselView
from particle_3d.particle82_provenance import atomic_json

p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args()
s=Scene(a.root);kind='01_full_network_rotating';frames=schedule(s,kind);index=40
before={pid:hashlib.sha256(array.tobytes()).hexdigest() for pid,array in s.arrays.items()}
images=[];receipts=[]
for _ in range(2):
    view=VesselView(s,wall_opacity=.12)
    im,receipt=compose(s,view,kind,index,len(frames),frames[index]);view.close()
    images.append(np.asarray(im));receipts.append(receipt)
result=dict(pixel_exact=bool(np.array_equal(*images)),frame_receipt_exact=receipts[0]==receipts[1],
    saved_state_unmodified=all(hashlib.sha256(s.arrays[pid].tobytes()).hexdigest()==sha for pid,sha in before.items()),
    input_arrays_readonly=all(not array.flags.writeable for array in s.arrays.values()),
    scope='FRESH_WINDOWS_IN_SAVED_SERVER_RENDER_ENVIRONMENT_NOT_A_CROSS_GPU_GUARANTEE',
    source_scene_sha256=s.sha256,pixel_sha256=[hashlib.sha256(im.tobytes()).hexdigest() for im in images])
result['all_pass']=all(result[k] for k in ['pixel_exact','frame_receipt_exact','saved_state_unmodified','input_arrays_readonly'])
atomic_json(Path(a.root)/'data/render_determinism.json',result)
if not result['all_pass']:raise ValueError('Saved replay determinism failed')
print(result)
