"""One deterministic, resolved three-lobe mural solid on the actual Y wall."""
import json,time
import numpy as np
import pyvista as pv
from scipy.ndimage import maximum_filter
from .config import ROOT,dump
from .y_bifurcation_sources import YSources,digest
from .y_mural_clot_geometry import CONFIG as BASE,build
from .y_surface_coordinates import validate_and_map

TASK='y_clot_realistic_morphology_wall_contact_and_flow_transport_v1'
OUT=ROOT/'verification'/TASK
RESULTS=ROOT/'results'/TASK
V2=ROOT/'verification/y_bifurcation_irregular_clot_and_initial_fracture_v1'
LIBDIR=ROOT/'tools/lammps-native-y-wall-shared-v1-build'
SOURCE=ROOT/'tools/lammps-contact-onset-v1'
DEPS=ROOT/'tools/native-y-wall-v1-python-deps'
CONFIG={**BASE,'geometry_id':'CURVED_MURAL_MULTILOBE_CLOT_V3','a_m':.00068,'b_m':.00062,'H_m':.00065,
        'geometry_seed':20261008,'randomness_used':False,'envelope_power':1.15,
        'lobe_centers_axial_arc_m':[[-.00028,-.00015],[0.,.00029],[.00031,-.00011]],
        'lobe_sigma_m':[.00014,.00013,.00014],'lobe_weights':[1.,.94,.86],
        'platform_height_m':.00044,'lobe_height_m':.00043,
        'dt_s':.5e-6,'duration_s':.0011605,'explicit_steps':4660,'dump_every':20,'mpi_ranks':2,
        'numerical_drag_rate_s_inv':500.,'material_id':'LPS_ELASTIC_MAPPING_FROM_LIU2024'}

def raw_profile(x,c):
    p=np.atleast_2d(x);a=p[:,0]-c['x0_m'];b=c['Rv_m']*np.arctan2(p[:,1],p[:,2]);u=a/c['a_m'];v=b/c['b_m']
    U=u-.10*u*u+.05;V=(v-.10*(1-u*u)-.055*u)/(1+.10*np.sin(np.pi*u));q2=U*U+V*V
    bump=sum(w*np.exp(-((a-cc[0])**2+(b-cc[1])**2)/(2*sig**2)) for cc,sig,w in zip(c['lobe_centers_axial_arc_m'],c['lobe_sigma_m'],c['lobe_weights']))
    platform=c['platform_height_m']*(1-.09*np.tanh(2*u))
    return q2,np.maximum(1-q2,0)**c['envelope_power']*(platform+c['lobe_height_m']*bump)

def profile(x,c):
    q,h=raw_profile(x,c);return q,h/c['height_normalization']

def freeze():
    assert not (OUT/'V3_GEOMETRY_CONFIG.json').exists()
    c=dict(CONFIG);a,b=np.meshgrid(np.linspace(-.00085,.00085,601),np.linspace(-.0008,.0008,601),indexing='ij')
    p=np.c_[c['x0_m']+a.ravel(),c['Rv_m']*np.sin(b.ravel()/c['Rv_m']),c['Rv_m']*np.cos(b.ravel()/c['Rv_m'])]
    h=raw_profile(p,c)[1].reshape(a.shape);c['height_normalization']=float(h.max()/c['H_m'])
    h/=c['height_normalization'];peaks=(h==maximum_filter(h,size=51))&(h>.00045)
    c['analytic_local_maxima']=[dict(axial_offset_m=float(x),circumferential_arc_m=float(y),height_m=float(z)) for x,y,z in zip(a[peaks],b[peaks],h[peaks])]
    dump(OUT/'V3_GEOMETRY_CONFIG.json',c)
    (OUT/'V3_MORPHOLOGY_DEFINITION.md').write_text('''# CURVED_MURAL_MULTILOBE_CLOT_V3
Single deterministic engineering morphology. Not segmented patient anatomy.
Coordinates: a=x-x0, b=Rv atan2(y,z), u=a/0.68mm, v=b/0.62mm.
U=u-.10u^2+.05; V=[v-.10(1-u^2)-.055u]/[1+.10sin(pi u)]; q^2=U^2+V^2.
Height raw=max(1-q^2,0)^1.15 * [0.44mm*(1-.09 tanh(2u)) + 0.43mm*sum(w_j exp(-((a-a_j)^2+(b-b_j)^2)/(2sigma_j^2)))].
Three actual lobe centers/weights/sigmas are in CONFIG; sigma130-140um, FWHM306-330um (6+ dx). A broad base preserves connectivity. Height and derivative vanish at the footprint edge. Thicker upstream platform and unequal lobe weights form an asymmetric shoulder and downstream taper.
One predetermined scalar height normalization sets max height0.65mm. No random draws. Native WALL triangles are clipped/subdivided; height goes inward from actual wall points. Fresh global50um cubic mask, 27 support points and generic degree<2 whisker removal; new typed bonds.
The analytic definition is frozen before assembly. At most one explicit pre-dynamics morphology correction is allowed if the normal-scale preview has no visible lobes; never change after observing fracture.
''')
    np.savez_compressed(OUT/'V3_ANALYTIC_HEIGHT_MAP.npz',axial_offset_m=a,circumferential_arc_m=b,height_m=h)
    dump(OUT/'V3_MORPHOLOGY_FREEZE.json',dict(time_epoch_s=time.time(),config_sha256=digest(OUT/'V3_GEOMETRY_CONFIG.json'),definition_sha256=digest(OUT/'V3_MORPHOLOGY_DEFINITION.md'),analytic_peak_count=int(peaks.sum()),dynamic_solves_started=0))
    print('FROZEN_ANALYTIC_MAXIMA',c['analytic_local_maxima'],flush=True)

def geometry():
    assert not (OUT/'GEOMETRY_SUMMARY.json').exists()
    c=json.loads((OUT/'V3_GEOMETRY_CONFIG.json').read_text());s=YSources();g,c,solid,top,attach=build(s,c,out=OUT,results=RESULTS,profile_fn=profile)
    fit,clear=validate_and_map(s,g,c,solid,top,attach,out=OUT)
    n=g['clot_count'];ids=np.flatnonzero(g['surface']);norm=g['normals'][ids];t=np.tile([1.,0.,0.],(len(ids),1))-norm[:,0,None]*norm
    z=dict(np.load(OUT/'SURFACE_NODES_AND_NORMALS.npz'));z.update(axial_tangent_xyz=t/np.linalg.norm(t,axis=1)[:,None],circumferential_tangent_xyz=np.cross(norm,t/np.linalg.norm(t,axis=1)[:,None]),axial_projected_xyz=t,semantics='V3_INITIAL_EXPOSED_SURFACE; CURRENT_POSITION_GAUSSIAN_ENVELOPE; REFERENCE_NORMAL_PROJECTION')
    np.savez_compressed(OUT/'V3_SURFACE_NODES_AND_NORMALS.npz',**z)
    np.savez_compressed(OUT/'V3_MATERIAL_POINTS_AND_BONDS.npz',**{k:v for k,v in g.items() if k!='info'})
    old=json.loads((V2/'GEOMETRY_SUMMARY.json').read_text());info=g['info']
    for k in list(info):
        if k.startswith('old_hemisphere_'):info.pop(k)
    info.update(exposed_nodes=len(ids),relative_volume_to_V2=info['continuous_mesh_volume_m3']/old['continuous_mesh_volume_m3']-1,discrete_relative_continuous_error=info['clot_volume_m3']/info['continuous_mesh_volume_m3']-1)
    dump(OUT/'GEOMETRY_SUMMARY.json',info);dump(OUT/'V2_V3_GEOMETRY_COMPARISON.json',dict(V2=old,V3=info))
    from .y_mural_clot_geometry import locator,hit
    loc=locator(s.exterior);bad=sum(hit(loc,g['X'][a],g['X'][b])[0] for a,b in g['bonds'][g['bond_types']==3]);assert not bad
    dump(OUT/'V3_WALL_FIT_AND_LUMEN_CHECK.json',dict(**fit,clearance=clear,exact_wall_chords_crossing=bad,exposed_attachment_areas_disjoint=True,geometry=info,analytic_lobes=c['analytic_local_maxima'],gate_A='GEOMETRY_CHECKS_PASS_PENDING_PREVIEW'))
    print('V3_GEOMETRY',info,flush=True)

def load_geometry():
    g=dict(np.load(OUT/'V3_MATERIAL_POINTS_AND_BONDS.npz'));g['clot_count']=int(g['clot_count']);g['info']=json.loads((OUT/'GEOMETRY_SUMMARY.json').read_text());return g

if __name__=='__main__':
    import sys
    {'freeze':freeze,'geometry':geometry}[sys.argv[1]]()
