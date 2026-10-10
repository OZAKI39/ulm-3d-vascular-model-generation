"""Four bounded native LAMMPS cases; zero-load gate; existing event/lineage analysis."""
import json, shutil, time, csv
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from .config import ROOT,dump
from .y_irregular_morphology import OUT,RESULTS,EXE,load_geometry
from .y_bifurcation_sources import digest
from .coherent_failure_card import failure_card,elastic_mapping
from .y_scene_task import run0 as initialize
from .calibrate_failure import prelude,dumps
from .current_position_provider import pulse
from .failure_map import analyze_case,csv_out
from .detached_fragment_metrics import analyze as lineage,shape_metrics
from . import runner

CASES=[('C0_ZERO_LOAD_DYNAMIC','LOCAL_GAUSSIAN_SHEAR',0.),('C1_GAUSSIAN_A100','LOCAL_GAUSSIAN_SHEAR',100.),('C2_GAUSSIAN_A200','LOCAL_GAUSSIAN_SHEAR',200.),('C3_OPPOSING_A100','OPPOSING_SHEAR',100.)]

def tangent(g):
    n=g['normals'];t=np.tile([1.,0.,0.],(len(n),1))-n[:,0,None]*n;t[~g['surface']]=0;return t

def shape(X,r,g):
    c=np.asarray(r['load_centers_m']);w=r['gaussian_sigma_m']
    a=np.exp(-np.sum((X-c[0])**2,axis=1)/(2*w*w))
    if r['pattern_id']=='OPPOSING_SHEAR':a-=np.exp(-np.sum((X-c[1])**2,axis=1)/(2*w*w))
    return a[:,None]*tangent(g)

def freeze():
    assert not (OUT/'CASE_CONFIG_FROZEN.json').exists()
    g=load_geometry();c=json.loads((OUT/'PLACEMENT_CONFIGURATION.json').read_text())
    prov=json.loads((ROOT/'verification/lammps_contact_transfer_and_fragment_size_v1/CONTACT_VARIANT_PROVENANCE.json').read_text())
    assert digest(EXE)==prov['executable_sha256']=='7cba88869fb61531b4e5fe5bbcf3b129f80cc1d9c19c0a4761d55e69b67cf95b'
    dump(OUT/'MATERIAL_CONTACT_IDENTITY.json',dict(material_card=failure_card(),elastic_mapping=elastic_mapping(),contact=prov,actual_parameters={k:c[k] for k in ['K_Pa','G_Pa','rho_kg_m3','s00_bulk','s00_interface','alpha','dx_m','delta_m','anchor_layers','contact_distance_factor','numerical_drag_rate_s_inv']},geometry_strength_validation='NOT_ESTABLISHED'))
    top=pv.read(RESULTS/'paraview/exposed_surface.vtp');pp=top.points
    uv=np.c_[(pp[:,0]-c['x0_m'])/c['a_m'],c['Rv_m']*np.arctan2(pp[:,1],pp[:,2])/c['b_m']]
    def point(u,v):return pp[np.argmin(np.sum((uv-[u,v])**2,axis=1))].tolist()
    common={**c,'load_source':'ENGINEERING_ANALYTIC_STREAMING_LOAD','load_provider':'CURRENT_POSITION_GAUSSIAN_WITH_V2_REFERENCE_TANGENT_PROJECTION',
            'gaussian_sigma_m':.00028,'T_cross_s':c['H_m']/np.sqrt(c['G_Pa']/c['rho_kg_m3']),
            'duration_s':.0011605,'dt_s':.5e-6,'explicit_steps':4660,'dump_every':20,'mpi_ranks':2,
            'map_stage':'Y_V2_INITIAL_EXPLORATORY_FRACTURE','nominal_duration_factor':None,
            'bulk_macro_strength_Pa':6780.9,'interface_macro_strength_Pa':1957.41,
            'macro_strength_scope':'HISTORICAL_FLAT_CALIBRATION_PROXY_NOT_CURVED_VALIDATION',
            'material_id':'LPS_ELASTIC_MAPPING_FROM_LIU2024','numerical_drag_rate_s_inv':500.,'EXPLORATORY_ONLY':True}
    common['duration_over_Tcross']=common['duration_s']/common['T_cross_s']
    requests=[{**common,'case_id':name,'pattern_id':pat,'tau_budget_Pa':amp,'load_centers_m':[point(-.22,.04)] if pat=='LOCAL_GAUSSIAN_SHEAR' else [point(-.38,.08),point(.38,.08)]} for name,pat,amp in CASES]
    dump(OUT/'CASE_CONFIG_FROZEN.json',dict(cases=requests,maximum_formal_dynamics=4,run0_count=1,absolute_pulse_not_new_2Tcross=True,initial_area_and_projection_frozen=True,force_contract='F_i=A_i P_i^V2 e_x g(x_i(t)) scale pulse(t); current-position spatial envelope; reference surface projection, not updated deformed-surface tangency',persistent_definition='Historical exact-membership >=3 consecutive saved frames; >=20us observed span here; not long-term biological persistence',freeze_time_epoch_s=time.time()))
    dump(OUT/'FROZEN_INPUT_HASHES.json',{f:digest(OUT/f) for f in ['GEOMETRY_V2_CONFIG.json','MORPHOLOGY_DEFINITION.md','CASE_CONFIG_FROZEN.json','Y_CLOT_V2_GEOMETRY.npz','Y_CLOT_V2_SURFACE_NODES_AND_NORMALS.npz','Y_CLOT_V2_BONDS_AND_ANCHORS.npz']})
    (OUT/'CLASSIFICATION_RULES.md').write_text('''Engineering categories in precedence order: NO_FAILURE only if bulk and interface broken fractions exactly zero. PULVERIZATION if >=50% initial clot nodes in 1-2 node components. WHOLE_CLOT_DETACHMENT if total attached mass zero and largest detached component >=80%. HYBRID if detached >=34-node components coexist with >=1% unresolved debris. DETACHED_RESOLVED_FRAGMENTS if any detached >=34-node component. ATTACHED_CORE_WITH_LOCAL_DEBRIS_SHEDDING if attached core and detached material coexist. Otherwise ATTACHED_DAMAGE_WITHOUT_DETACHMENT or UNCERTAIN. FI_free retains historical 1-largest clot-only component/N (includes attached core in maximum); >=3-node count is historical engineering graph count, NOT a large-fragment label. Large threshold34/114 nodes =200/300um volume-equivalent diameter. Exact-membership lineage resets on every split; persistent >=3 saved frames is an engineering convention only.\n''')

def run0():
    initialize(out=OUT,geometry_loader=load_geometry)
    assert json.loads((OUT/'LAMMPS_RUN0_CHECK.json').read_text())['maximum_initial_force_N']==0
    shutil.copy2(OUT/'Y_CLOT_LAMMPS.data',OUT/'Y_CLOT_V2_LAMMPS.data')

def prepare(r,g):
    p=OUT/'cases'/r['case_id'];p.mkdir(parents=True,exist_ok=False)
    shutil.copy2(OUT/'Y_CLOT_V2_LAMMPS.data',p/'body.data')
    sh=shape(g['X'],r,g);area=float(g['areas'].sum());factor=r['tau_budget_Pa']*area/np.sum(np.linalg.norm(sh,axis=1)*g['areas']);F=sh*g['areas'][:,None]*factor
    load=dict(scale_Pa=factor,reference_peak_traction_Pa=float(np.linalg.norm(sh,axis=1).max()*factor),reference_surface_area_m2=area,initial_integrated_force_budget_N=float(np.linalg.norm(F,axis=1).sum()),reference_force_N=F,reference_net_force_N=F.sum(0),reference_moment_about_footpoint_Nm=np.cross(g['X']-r['actual_wall_footpoint_m'],F).sum(0),normalization='Frozen reference mean magnitude budget, no dynamic renormalization',direction='New V2 initial exposed normals projected axial basis; frozen in global frame',field='Spatial Gaussian envelope sampled at current xyz every step; no particle-following center',new_fracture_surfaces_loaded=False)
    dump(p/'request.json',r);dump(p/'LOAD_PROVIDER.json',load)
    ids=np.flatnonzero(g['surface']);tt=tangent(g)
    text=prelude(g,r,r['dump_every']).replace('official unmodified stable_30Sep2026','verified existing modified 1.0dx stable_30Sep2026').replace('pair_style bpm/peri\n','pair_style bpm/peri contact_distance_factor 1.0\n')
    text+='fix tangent_props all property/atom d_tx d_ty d_tz ghost yes\nset group all d_tx 0\nset group all d_ty 0\nset group all d_tz 0\ngroup surface id '+' '.join(map(str,ids+1))+'\n'
    for i in ids:text+=f'set atom {i+1} d_area {g["areas"][i]:.17g} d_tx {tt[i,0]:.17g} d_ty {tt[i,1]:.17g} d_tz {tt[i,2]:.17g}\n'
    w=r['gaussian_sigma_m']
    for j,center in enumerate(r['load_centers_m']):
        exponent='+'.join(f'({a}-({c:.17g}))^2' for a,c in zip('xyz',center))
        text+=f'variable env{j} atom exp(-({exponent})/{2*w*w:.17g})\n'
    text+='variable spatial atom '+('v_env0-v_env1' if r['pattern_id']=='OPPOSING_SHEAR' else 'v_env0')+'\n'
    T=r['duration_s'];text+=f'variable amp equal ternary(time<{T/4:.17g},0.5*(1-cos(4*PI*time/{T:.17g})),ternary(time<{.75*T:.17g},1.0,ternary(time<{T:.17g},0.5*(1+cos(4*PI*(time/{T:.17g}-0.75))),0.0)))\n'
    for a in 'xyz':text+=f'variable load_f{a} atom d_area*{factor:.17g}*v_amp*v_spatial*d_t{a}\n'
    text+=f'''fix stream surface addforce v_load_fx v_load_fy v_load_fz
fix drag clot viscous {r['rho_kg_m3']*r['dx_m']**3*r['numerical_drag_rate_s_inv']:.17g}
fix fixed anchors setforce 0 0 0
fix integrate clot nve
timestep {r['dt_s']:.17g}
variable force_mag atom sqrt(v_load_fx^2+v_load_fy^2+v_load_fz^2)
variable local_tau atom v_force_mag/(d_area+1e-100)
compute force_budget surface reduce sum v_force_mag v_load_fx v_load_fy v_load_fz
compute peak_tau surface reduce max v_local_tau
variable actual_budget equal c_force_budget[1]/{area:.17g}
fix force_history all ave/time 1 1 1 v_actual_budget c_peak_tau c_force_budget[1] c_force_budget[2] c_force_budget[3] c_force_budget[4] file force_history.dat
thermo 100
thermo_style custom step time atoms bonds v_amp v_actual_budget ke
thermo_modify lost error lost/bond error flush yes format float %.15g
timer full
'''+dumps(r['dump_every'],r['explicit_steps'],'v_load_fx v_load_fy v_load_fz')
    (p/'in.lammps').write_text(text)
    return dict(path=p,request=r,geometry=g,load=load,steps=r['explicit_steps'])

def classify(m):
    if m['bulk_broken_fraction']==0 and m['interface_broken_fraction']==0:return 'NO_FAILURE'
    if m['debris_fraction']>=.5:return 'PULVERIZATION'
    if m['total_attached_fraction']==0 and m['largest_detached_nodes']>=.8*CURRENT_N:return 'WHOLE_CLOT_DETACHMENT'
    if m['largest_detached_nodes']>=34:return 'HYBRID' if m['debris_fraction']>=.01 else 'DETACHED_RESOLVED_FRAGMENTS'
    if 0<m['total_attached_fraction']<1:return 'ATTACHED_CORE_WITH_LOCAL_DEBRIS_SHEDDING'
    if m['total_attached_fraction']==1:return 'ATTACHED_DAMAGE_WITHOUT_DETACHMENT'
    return 'UNCERTAIN'

def analyze(case):
    global CURRENT_N
    started=time.perf_counter();g=case['geometry'];CURRENT_N=g['clot_count'];p=case['path'];r=case['request'];n=CURRENT_N
    result=analyze_case(case,shape_function=lambda X,r:shape(X,r,g))
    # Native per-timestep integrated external force; also saved-frame moment integral.
    fh=np.loadtxt(p/'force_history.dat');times=np.r_[0,fh[:,0]*r['dt_s']];vals=np.vstack([np.zeros(fh.shape[1]-1),fh[:,1:]])
    result.update(actual_peak_local_traction_Pa=float(vals[:,1].max()),peak_actual_mean_traction_budget_Pa=float(vals[:,0].max()),integral_sum_force_magnitudes_Ns=float(np.trapezoid(vals[:,2],times)),net_impulse_Ns=np.trapezoid(vals[:,3:6],times,axis=0),force_history_resolution_s=r['dt_s'])
    dump(p/'result.json',result)
    summ=lineage(p,p/'lineage',r['case_id'],g,1.,shape_min_nodes=34,classify_fn=classify)
    tr=np.load(p/'accepted_trajectory.npz');rows=list(csv.DictReader((p/'lineage/timeseries.csv').open()));oldrows=list(csv.DictReader((p/'fragment_metrics.csv').open()))
    moments=[];current_clearance=[]
    from .y_bifurcation_sources import YSources
    from .y_mural_clot_geometry import distances
    src=YSources()
    for k,row in enumerate(rows):
        labels=tr['component_id'][k];att=tr['attached'][k];count=np.unique(labels[~att],return_counts=True)[1]
        row.update(FI_free=oldrows[k]['FI_free'],historical_ge3_component_count=oldrows[k]['resolved_fragment_count'],detached_component_count=len(count),detached_mass_fraction=float((~att).mean()),detached_1_2node_debris_fraction=float(count[count<=2].sum()/n),actual_mean_traction_budget_Pa=oldrows[k]['actual_mean_traction_budget_Pa'])
        X=tr['position_m'][k];F=shape(X,r,g)*g['areas'][:,None]*case['load']['scale_Pa']*pulse(int(tr['step'][k])*r['dt_s'],r['duration_s']);moments.append(np.cross(X-r['actual_wall_footpoint_m'],F).sum(0))
        outside=distances(src.exterior,X[:n])>2e-9;row['material_nodes_outside_visual_vessel']=int(outside.sum());current_clearance.append(int(outside.sum()))
    csv_out(p/'lineage/timeseries.csv',rows)
    comp=list(csv.DictReader((p/'lineage/components.csv').open()));patch=pv.read(RESULTS/'paraview/attachment_surface.vtp');tree=cKDTree(patch.points)
    for row in comp:
        row['COM_distance_to_original_attachment_patch_vertices_m']=float(tree.query([float(row['COM_'+a+'_m']) for a in 'xyz'])[0])
    csv_out(p/'lineage/components.csv',comp)
    ev=np.load(p/'broken_bonds.npz');first={}
    for t,label in [(1,'bulk'),(2,'interface')]:
        e=ev['events'][g['bond_types'][ev['initial_bond_indices']]==t]
        first[label]=None if not len(e) else dict(time_s=float(e[0,1]*r['dt_s']),node_ids=e[0,2:4].astype(int),position_m=e[0,4:7],reference_position_m=e[0,7:10])
    result.update(final_classification=summ['final']['failure_mode'],first_break=first,lineage_summary=summ,maximum_nodes_outside_visual_vessel=max(current_clearance),angular_impulse_about_footpoint_Nms=np.trapezoid(moments,tr['step']*r['dt_s'],axis=0),angular_impulse_sampling_s=r['dt_s']*r['dump_every'],postprocess_total_wall_s=time.perf_counter()-started)
    dump(p/'result.json',result);print('CASE_COMPLETE',r['case_id'],result['final_classification'],'bulk',result['bulk_broken_fraction'],'interface',result['interface_broken_fraction'],'attached',result['attached_clot_volume_fraction'],flush=True)
    return result

def run(case_index):
    frozen=json.loads((OUT/'CASE_CONFIG_FROZEN.json').read_text());assert 0<=case_index<4
    for name,sha in json.loads((OUT/'FROZEN_INPUT_HASHES.json').read_text()).items():assert digest(OUT/name)==sha
    assert json.loads((OUT/'LAMMPS_RUN0_CHECK.json').read_text())['status']=='PASS'
    if case_index:assert json.loads((OUT/'ZERO_LOAD_GATE.json').read_text())['status']=='PASS'
    assert digest(EXE)==json.loads((OUT/'MATERIAL_CONTACT_IDENTITY.json').read_text())['contact']['executable_sha256']
    g=load_geometry();case=prepare(frozen['cases'][case_index],g);runner.LMP=EXE
    perf=runner.run_lammps(case);perf['scope']='Actual Y V2 mural clot case';dump(case['path']/'PERFORMANCE.json',perf)
    result=analyze(case)
    if case_index==0:
        a=np.load(case['path']/'accepted_trajectory.npz');dis=float(abs(a['position_m']-g['X']).max());speed=float(abs(a['velocity_m_s']).max())
        okay=result['bulk_broken_fraction']==0 and result['interface_broken_fraction']==0 and dis<1e-12 and speed<1e-10 and result['actual_peak_local_traction_Pa']==0
        dump(OUT/'ZERO_LOAD_GATE.json',dict(status='PASS' if okay else 'FAIL',max_displacement_m=dis,max_speed_m_s=speed,no_broken_bonds=result['bulk_broken_fraction']==result['interface_broken_fraction']==0,finite_states=result['validation']['finite_states'],steps=4660))
        assert okay,'Zero-load gate failed; stop loaded cases'

if __name__=='__main__':
    import sys
    if sys.argv[1]=='freeze':freeze()
    elif sys.argv[1]=='run0':run0()
    elif sys.argv[1]=='run':run(int(sys.argv[2]))
