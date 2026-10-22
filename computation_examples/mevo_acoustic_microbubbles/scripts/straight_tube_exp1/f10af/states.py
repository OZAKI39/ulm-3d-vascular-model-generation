"""Frozen random states; synthetic pattern tests never enter this library."""
import numpy as np,h5py
from scipy.spatial.distance import pdist
from .common import *
from .field import response_context
def dispersed(n,seed,scale=1.):
    rng=np.random.default_rng(seed);x=[]
    while len(x)<n:
        q=rng.uniform([-.0005,-.0004,-.00002],[.0005,.0004,.00002])*scale
        if not x or np.linalg.norm(np.asarray(x)-q,axis=1).min()>8*R:x.append(q)
    return np.asarray(x)
def main():
    verify();c=response_context();states={};cfg=read(CONFIG)
    for n in [32,64,128,200]:states[f'N{n}_wide_r0']=dispersed(n,cfg['seed_discovery']+n)
    for j,seed in enumerate(cfg['seeds_validation']):states[f'N64_wide_v{j}']=dispersed(64,seed)
    states['N64_loose_r0']=dispersed(64,462030,.2)
    with h5py.File(AE/'INITIAL_STATES.h5') as h:states['N64_compact_reference']=h['S2_r0'][:]
    q=dispersed(32,462040);q[:2]=[[-(R+.4*R),0,0],[(R+.4*R),0,0]];states['N32_contact_transition']=q
    q=dispersed(32,462041);q[:,1]*=.2
    q[:,2]=np.sqrt((.0008-R-5e-6)**2-q[:,1]**2);states['N32_prescribed_top']=q
    index=[];bounds=np.array(cfg['field']['bounds_m']);margin=2*cfg['field']['coarse_step_m']
    with h5py.File(BASE/'INITIAL_STATES.h5','x') as h:
        for name,x in states.items():
            checks=[c['cavity'].query(z,R) for z in x];cache=bool(np.all(x>bounds[0]+margin) and np.all(x<bounds[1]-margin));gap=float(pdist(x).min()-2*R)
            record=dict(id=name,N=len(x),geometry_valid=all(q['valid'] for q in checks) and gap>=0,cache_valid=cache,minimum_pair_gap_m=gap,minimum_wall_gap_m=min(q['finite_radius_gap_m'] for q in checks),extent_m=np.ptp(x,axis=0),height_m=[x[:,2].min(),x[:,2].max()],source='F10AE connected demonstration' if 'compact' in name else 'Frozen random rejection sampling, not acoustically formed',scope='PRESCRIBED_NEAR_TOP_WALL_INITIAL_CONDITION; dynamic cache NOT_BUILT' if 'top' in name else 'SIMULATION_3D_TRUTH')
            assert record['geometry_valid'];h[name]=x;index.append(record)
    write(BASE/'INITIAL_STATE_INDEX.json',dict(states=index,sha256=sha(BASE/'INITIAL_STATES.h5'),seeds=cfg['seeds_validation'],identity=identity()));print('STATES',len(index),flush=True)
if __name__=='__main__':main()
