from pathlib import Path
import sys,json
from types import SimpleNamespace
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'particle_3d/src'));sys.path.insert(0,str(ROOT/'particle_3d/scripts'))
from prepare_rbc_mb_coflow import OUT,RADIUS,MEAN,DURATION,sha
from particle_3d.normal_rbc_encounter import poiseuille,MB_RADIUS,RBC_RADIUS
from particle_3d.microbubble import MicrobubbleState
from particle_3d.integrator import advance_single_microbubble
from particle_3d.rbc import RBCGeometry,RBCState
from particle_3d.rbc_distribution import sample_rbc_geometries
from particle_3d.rbc_integrator import advance_single_rbc

def load():
    return json.loads((OUT/'data/SCENE.json').read_text()),np.load(OUT/'data/coflow.npz')

class AnalyticalField:
    def sample(self,p):
        gradient=np.zeros((3,3));gradient[0,1:]=-4*MEAN*np.asarray(p)[1:]/RADIUS**2
        curl=np.array([0.,gradient[0,2],-gradient[0,1]])
        return SimpleNamespace(inside_lumen=True,tetra_id=0,velocity_m_s=poiseuille(p,RADIUS,MEAN),
            velocity_gradient_s_inv=gradient,vorticity_s_inv=curl)

def test_closed_form_transport_matches_existing_particle_integrators():
    scene,data=load();field=AnalyticalField();g=RBCGeometry.from_population(sample_rbc_geometries(41311,2026092002),41310)
    for j in [0,3,8,9,20,36]:
        r=scene['records'][j];p=np.array(r['initial_m']);z=np.zeros(3)
        if r['species']=='RBC':
            state=RBCState(r['id'],p,[1,0,0,0],z,z,g);end=advance_single_rbc(state,field,DURATION)
            np.testing.assert_allclose(end.short_axis,[0,0,1],rtol=0,atol=1e-15)
        else:
            state=MicrobubbleState(r['id'],p,r['radius_m'],z,z);end=advance_single_microbubble(state,field,DURATION)
        np.testing.assert_allclose(end.position_m,data['positions_m'][-1,j],rtol=0,atol=1e-18)

def test_shared_clock_real_sizes_and_complete_unmodified_arrays():
    scene,data=load();t=data['time_s'];x=data['positions_m']
    assert x.shape==(432,37,3) and scene['counts']=={'RBC':9,'MB':28}
    assert t[0]==0 and t[-1]==DURATION and np.all(np.diff(t)>0)
    assert sha(OUT/'data/coflow.npz')==scene['array_sha256']
    assert sha(OUT/'data/normal_rbc_si.vtp')==scene['shape_source_sha256']
    for j,r in enumerate(scene['records']):
        np.testing.assert_array_equal(x[:,j,1:],np.tile(r['initial_m'][1:],(432,1)))
        assert r['radius_m']==(RBC_RADIUS if r['species']=='RBC' else MB_RADIUS)
    assert not scene['pair_contact'] and not scene['hydrodynamic_pair_coupling']

def test_whole_interval_wall_and_pair_clearance_without_repulsion():
    scene,data=load();records=scene['records'];x0=data['positions_m'][0];v=data['velocity_m_s']
    z_half=scene['original_geometry']['maximum_thickness_um']*.5e-6
    for j,r in enumerate(records):
        # Fixed normal disc: actual meridional extent, bounded by a larger box.
        radial=np.hypot(abs(x0[j,1])+RBC_RADIUS,z_half) if r['species']=='RBC' else np.linalg.norm(x0[j,1:])+MB_RADIUS
        assert radial<RADIUS
    for i,a in enumerate(records):
        for j in range(i):
            b=records[j];d=x0[i]-x0[j];w=v[i]-v[j]
            t=0. if w@w==0 else np.clip(-d@w/(w@w),0,DURATION)
            nearest=d+t*w
            if a['species']==b['species']:
                assert np.linalg.norm(nearest)>a['radius_m']+b['radius_m']
            else:
                # Separate bounding slabs certify finite-size clearance at ALL x.
                assert abs(d[1])>RBC_RADIUS+MB_RADIUS or abs(d[2])>z_half+MB_RADIUS

def test_observation_window_contains_both_species_at_every_shared_time():
    scene,data=load()
    for pos in data['positions_m']:
        counts={'RBC':0,'MB':0}
        for p,r in zip(pos,scene['records']):
            if p[0]+r['radius_m']>=-20e-6 and p[0]-r['radius_m']<=20e-6:counts[r['species']]+=1
        assert counts['RBC']>=2 and counts['MB']>=4
