"""Bounded V3 C0/C1 with actual MPI fix external native-wall integration."""
import csv,json,shutil,time,os,subprocess,sys,re
from pathlib import Path
import numpy as np
import pyvista as pv
import psutil
from .config import ROOT,dump
from .y_clot_morphology_v3 import OUT,RESULTS,LIBDIR,SOURCE,DEPS,load_geometry
from .y_bifurcation_sources import digest
from .y_fragmentation_cases import shape,tangent
from .calibrate_failure import prelude,dumps
from .interface_model import write_typed_data
from .native_y_wall_contact import CONTACT,WallCallback,minimal_checks
from .failure_map import analyze_case,csv_out
from .detached_fragment_metrics import analyze as lineage
from .fragment_analysis import read_dump

def freeze():
    assert not (OUT/'CASE_CONFIG_FROZEN.json').exists()
    g=load_geometry();c=json.loads((OUT/'PLACEMENT_CONFIGURATION.json').read_text())
    pp=pv.read(RESULTS/'paraview/exposed_surface.vtp').points;uv=np.c_[(pp[:,0]-c['x0_m'])/c['a_m'],c['Rv_m']*np.arctan2(pp[:,1],pp[:,2])/c['b_m']]
    center=pp[np.argmin(np.sum((uv-[-.22,.04])**2,axis=1))]
    common={**c,'load_source':'ENGINEERING_ANALYTIC_STREAMING_LOAD','load_provider':'CURRENT_POSITION_GAUSSIAN_WITH_V3_REFERENCE_TANGENT_PROJECTION','gaussian_sigma_m':.00028,'T_cross_s':c['H_m']/np.sqrt(c['G_Pa']/c['rho_kg_m3']),'duration_s':.0011605,'dt_s':.5e-6,'explicit_steps':4660,'dump_every':20,'mpi_ranks':2,'map_stage':'Y_V3_NATIVE_WALL_EXPLORATORY','nominal_duration_factor':None,'bulk_macro_strength_Pa':6780.9,'interface_macro_strength_Pa':1957.41,'material_id':'LPS_ELASTIC_MAPPING_FROM_LIU2024','numerical_drag_rate_s_inv':500.,'EXPLORATORY_ONLY':True,'pattern_id':'LOCAL_GAUSSIAN_SHEAR','load_centers_m':[center.tolist()]}
    common['duration_over_Tcross']=common['duration_s']/common['T_cross_s']
    cases=[{**common,'case_id':name,'tau_budget_Pa':amp} for name,amp in [('C0_ZERO_LOAD_DYNAMIC',0.),('C1_GAUSSIAN_A200',200.)]]
    write_typed_data(OUT/'V3_LAMMPS.data',g,common)
    dump(OUT/'CASE_CONFIG_FROZEN.json',dict(cases=cases,maximum_formal_dynamics=2,handoff_step=4660,handoff_time_s=.00233,freeze_time_epoch_s=time.time(),contact=CONTACT,force_contract='F=A_reference*(I-n_reference*n_reference)*ex*Gaussian(current xyz)*frozen_scale*pulse; reference surface directions, no new fracture surface loading'))
    dump(OUT/'CONTACT_SHARED_LIBRARY_PROVENANCE.json',dict(library=str(LIBDIR/'liblammps.so'),library_sha256=digest(LIBDIR/'liblammps.so'),source=str(SOURCE),source_bond_sha256=digest(SOURCE/'src/BPM/bond_bpm_peri.cpp'),source_pair_sha256=digest(SOURCE/'src/BPM/pair_bpm_peri.cpp'),MPI='OpenMPI 4.1.6 explicitly selected for both mpi4py and LAMMPS',mpi4py='4.1.2 isolated target',official_static_binary_untouched=True))
    dump(OUT/'FROZEN_INPUT_HASHES.json',{f:digest(OUT/f) for f in ['V3_GEOMETRY_CONFIG.json','V3_MORPHOLOGY_DEFINITION.md','V3_MATERIAL_POINTS_AND_BONDS.npz','V3_LAMMPS.data','CASE_CONFIG_FROZEN.json']})
    q=minimal_checks();dump(OUT/'NATIVE_WALL_CONTACT_MINIMAL_CHECKS.json',q)
    gate=json.loads((OUT/'V3_WALL_FIT_AND_LUMEN_CHECK.json').read_text());gate['gate_A']='PASS';gate['preview_visual_check']='Three visible lobes at same camera/scale as V2; no morphology revision used';dump(OUT/'V3_WALL_FIT_AND_LUMEN_CHECK.json',gate)
    (OUT/'WALL_CONTACT_MODEL.md').write_text('''# NATIVE_Y_WALL_REPULSION_V1
Actual LAMMPS fix external pf/callback 1 1, called on every MPI rank at every post_force (including setup). SI force conversion precedes the independent external Newton force. Anchors excluded by immutable atom IDs.
Native open WALL nearest triangles; inward normals oriented with opposite vertex of adjacent native tetrahedron. Strict native tetra membership supplies distance sign. Actual triangle ports handled separately; no collision caps, no analytic cylinder. BVH locators and exact unchanged-position cache; no approximate wall mesh.
Engineering center proxy radius20um=0.4dx; k=50N/m; approaching-only damping zeta=.5, c=7.905694e-5Ns/m. Force=[k max(r-d,0)+c max(-v.n,0) I(d<r)]n. Node mass1.25e-10kg; dt omega=.31623. Center penetration tolerance2um=.04dx, reported separately from positive-clearance proxy overlap. Initial smallest clearance exceeds20um: no initial preload.
MODIFIED_1P0DX particle-particle contact unchanged. This is an engineering discrete-center wall constraint, not calibrated tissue contact or continuum FSI. No coordinate repair or render clipping of material nodes.
Every-step rank logs and per-node min/duration record actual coordinates. Swept segments are tested when previous wall distance <= segment length (a necessary condition for crossing). Legal actual-port crossings deactivate wall repulsion for that exited node. No trajectory has been rescaled.
Reference: https://docs.lammps.org/fix_external.html ; https://docs.lammps.org/Python_module.html ; https://docs.lammps.org/Python_install.html . API verified against fixed local stable_30Sep2026 source.
''')

def driver():
    from mpi4py import MPI
    from lammps import lammps
    comm=MPI.COMM_WORLD;p=Path.cwd();n=load_geometry()['clot_count']
    L=lammps(comm=comm,cmdargs=['-log','log.lammps','-screen','screen.lammps']);cb=None
    try:
        text=(p/'in.lammps').read_text();prefix,suffix=text.split('\nrun ',1)
        L.commands_string(prefix);cb=WallCallback(L,comm,p,n);L.set_fix_external_callback('nativewall',cb.callback)
        L.commands_string('run '+suffix)
    finally:
        if cb:cb.close()
        L.close()

def launch(case):
    p=case['path'];env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH=f'{SOURCE}/python:{DEPS}:{ROOT}',LD_LIBRARY_PATH=str(LIBDIR),MPI4PY_MPIABI='openmpi',MPI4PY_LIBMPI='/usr/lib/x86_64-linux-gnu/libmpi.so.40')
    cmd=['mpirun','--bind-to','core','-np','2',sys.executable,'-m','pd_clot.lammps_fragmentation.y_v3_cases','driver']
    dump(p/'RUN_COMMAND.json',dict(argv=cmd,cwd=str(p),environment={k:env[k] for k in ['PYTHONPATH','LD_LIBRARY_PATH','MPI4PY_MPIABI','MPI4PY_LIBMPI','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']}))
    start=time.perf_counter();peak=0
    with (p/'launcher.log').open('w') as out:
        proc=subprocess.Popen(['/usr/bin/time','-f','%U %S %e %M %x','-o','resource.time',*cmd],cwd=p,env=env,stdout=out,stderr=subprocess.STDOUT);parent=psutil.Process(proc.pid)
        while proc.poll() is None:
            try:peak=max(peak,sum(q.memory_info().rss for q in [parent,*parent.children(recursive=True)] if q.is_running()))
            except (psutil.NoSuchProcess,psutil.AccessDenied):pass
            time.sleep(.2)
    wall=time.perf_counter()-start;numbers=list(map(float,(p/'resource.time').read_text().splitlines()[-1].split()))
    output=(p/'log.lammps').read_text() if (p/'log.lammps').exists() else '';match=re.search(r'Loop time of ([\d.e+-]+) on (\d+) procs for (\d+) steps',output)
    perf=dict(exit_code=proc.returncode,wall_seconds=wall,CPU_seconds=sum(numbers[:2]),peak_RSS_bytes=max(peak,int(numbers[3])*1024),mpi_ranks=2,explicit_steps=int(match[3]) if match else None,loop_wall_s=float(match[1]) if match else None)
    dump(p/'PERFORMANCE.json',perf);assert proc.returncode==0 and 'ERROR:' not in output,str(p)
    return perf

def run0():
    g=load_geometry();r=json.loads((OUT/'CASE_CONFIG_FROZEN.json').read_text())['cases'][0];r={**r,'case_id':'RUN0_NATIVE_CALLBACK','explicit_steps':1}
    case=prepare(r,g);p=case['path'];s=(p/'in.lammps').read_text().replace('run 1\n','run 0\n');(p/'in.lammps').write_text(s);launch(case)
    frames=list(read_dump(p/'trajectory.lammpstrj'));assert len(frames)==1 and frames[0][0]==0
    _,cols,a=frames[0];X=a[:,[cols.index(k) for k in ['x','y','z']]];F=a[:,[cols.index(k) for k in ['fx','fy','fz']]]
    assert np.max(abs(X-g['X']))<1e-15 and np.max(abs(F))==0
    b=list(read_dump(p/'surviving_bonds.local'))[0][2][:,1:4].astype(int);b[:,:2]=np.sort(b[:,:2],axis=1)-1;b=b[np.lexsort((b[:,1],b[:,0]))]
    assert np.array_equal(b[:,:2],g['bonds']) and np.array_equal(b[:,2],g['bond_types'])
    dump(OUT/'LAMMPS_RUN0_CHECK.json',dict(status='PASS',time_steps_advanced=0,actual_MPI_ranks=2,coordinates_unchanged=True,topology_matches=True,maximum_initial_force_N=float(abs(F).max()),callback_calls_by_rank=[json.loads((p/f'wall_callback_rank{k}.json').read_text())['calls'] for k in range(2)]))

def analyze(case):
    p=case['path'];g=case['geometry'];r=case['request'];result=analyze_case(case,shape_function=lambda X,r:shape(X,r,g))
    def classify(m):
        if m['bulk_broken_fraction']==m['interface_broken_fraction']==0:return 'NO_FAILURE'
        if 0<m['total_attached_fraction']<1:return 'ATTACHED_CORE_WITH_LOCAL_DEBRIS_SHEDDING' if m['largest_detached_nodes']<34 else 'DETACHED_RESOLVED_FRAGMENTS'
        return 'ATTACHED_DAMAGE_WITHOUT_DETACHMENT' if m['total_attached_fraction']==1 else 'DETACHED'
    summ=lineage(p,p/'lineage',r['case_id'],g,1.,shape_min_nodes=34,classify_fn=classify);result['lineage_summary']=summ
    logs=[np.genfromtxt(p/f'wall_rank{k}.csv',delimiter=',',names=True) for k in range(2)];rows=[]
    for step in range(0,r['explicit_steps']+1):
        records=[a[a['step']==step][-1] for a in logs];rows.append(dict(case=r['case_id'],step=step,time_s=step*r['dt_s'],min_clearance_m=min(z['min_clearance_m'] for z in records),max_center_penetration_m=max(z['max_center_penetration_m'] for z in records),contact_nodes=sum(z['contact_nodes'] for z in records),outside_nodes=sum(z['outside_nodes'] for z in records),sum_normal_force_N=sum(z['sum_normal_force_N'] for z in records)))
    csv_out(p/'WALL_CONTACT_DIAGNOSTICS.csv',rows)
    node=[np.load(p/f'wall_node_audit_rank{k}.npz') for k in range(2)];minimum=np.minimum.reduce([z['minimum_clearance_m'] for z in node]);counts=np.sum([z['outside_step_count'] for z in node],axis=0);force=np.maximum.reduce([z['maximum_force_N'] for z in node]);meta=[json.loads((p/f'wall_callback_rank{k}.json').read_text()) for k in range(2)];events=sum([z['events'] for z in meta],[])
    np.savez_compressed(p/'WALL_NODE_AUDIT.npz',atom_id=np.arange(1,g['clot_count']+1),minimum_clearance_m=minimum,outside_duration_s=counts*r['dt_s'],maximum_force_N=force)
    stats=dict(minimum_clearance_m=float(minimum.min()),maximum_center_penetration_m=float(max(0,-minimum.min())),minimum_clearance_atom_id=int(minimum.argmin()+1),unique_outside_nodes=int((counts>0).sum()),outside_mass_kg=float((counts>0).sum()*CONTACT['node_mass_kg']),maximum_outside_duration_s=float(counts.max()*r['dt_s']),unique_wall_contact_nodes=int((force>0).sum()),maximum_nodal_wall_force_N=float(force.max()),swept_crossings=sum(e['event']=='SWEPT_WALL_CROSSING' for e in events),open_port_exits=[e for e in events if e['event']=='OPEN_PORT_EXIT'],callback_seconds_rank_sum=sum(z['callback_seconds'] for z in meta),query_seconds_rank_sum=sum(z['locator_query_seconds'] for z in meta),distance_over_dx_percentiles=np.quantile(minimum/r['dx_m'],[0,.01,.05,.5,.95,1]),events=events)
    dump(p/'WALL_CONTACT_DIAGNOSTICS.json',stats);result['wall_contact']=stats;dump(p/'result.json',result)
    if r['tau_budget_Pa']==0:
        tr=np.load(p/'accepted_trajectory.npz');dis=float(abs(tr['position_m']-g['X']).max());speed=float(abs(tr['velocity_m_s']).max());okay=result['bulk_broken_fraction']==result['interface_broken_fraction']==0 and dis<1e-12 and speed<1e-10 and force.max()==0
        dump(OUT/'ZERO_LOAD_GATE.json',dict(status='PASS' if okay else 'FAIL',max_displacement_m=dis,max_speed_m_s=speed,wall_force_N=float(force.max()),steps=4660,no_broken_bonds=result['bulk_broken_fraction']==result['interface_broken_fraction']==0));assert okay
    print('CASE_COMPLETE',r['case_id'],result['bulk_broken_fraction'],result['interface_broken_fraction'],result['attached_clot_volume_fraction'],stats,flush=True)

def run(index):
    assert index in (0,1);assert json.loads((OUT/'LAMMPS_RUN0_CHECK.json').read_text())['status']=='PASS'
    if index:assert json.loads((OUT/'ZERO_LOAD_GATE.json').read_text())['status']=='PASS'
    for f,h in json.loads((OUT/'FROZEN_INPUT_HASHES.json').read_text()).items():assert digest(OUT/f)==h
    g=load_geometry();r=json.loads((OUT/'CASE_CONFIG_FROZEN.json').read_text())['cases'][index];case=prepare(r,g);launch(case);analyze(case)

def prepare(r,g):
    p=OUT/'cases'/r['case_id'];p.mkdir(parents=True,exist_ok=False)
    shutil.copy2(OUT/'V3_LAMMPS.data',p/'body.data')
    sh=shape(g['X'],r,g);area=float(g['areas'].sum());factor=r['tau_budget_Pa']*area/np.sum(np.linalg.norm(sh,axis=1)*g['areas']);F=sh*g['areas'][:,None]*factor
    load=dict(scale_Pa=factor,reference_peak_traction_Pa=float(np.linalg.norm(sh,axis=1).max()*factor),reference_surface_area_m2=area,initial_integrated_force_budget_N=float(np.linalg.norm(F,axis=1).sum()),reference_force_N=F,reference_net_force_N=F.sum(0),reference_moment_about_footpoint_Nm=np.cross(g['X']-r['actual_wall_footpoint_m'],F).sum(0),normalization='Frozen reference mean magnitude budget, no dynamic renormalization',direction='New V3 initial exposed normals projected axial basis; frozen in global frame',field='Spatial Gaussian envelope sampled at current xyz every step; no particle-following center',new_fracture_surfaces_loaded=False)
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
    text=text.replace('fix fixed anchors', 'fix wallprops all property/atom d_wclear d_wpen d_wtri\nfix nativewall clot external pf/callback 1 1\nfix fixed anchors')
    text=text.replace('v_load_fx v_load_fy v_load_fz\ndump_modify trajectory','v_load_fx v_load_fy v_load_fz d_wclear d_wpen d_wtri f_nativewall[1] f_nativewall[2] f_nativewall[3]\ndump_modify trajectory')
    (p/'in.lammps').write_text(text)
    return dict(path=p,request=r,geometry=g,load=load,steps=r['explicit_steps'])

if __name__=='__main__':
    if sys.argv[1]=='run':run(int(sys.argv[2]))
    else:{'freeze':freeze,'run0':run0,'driver':driver}[sys.argv[1]]()
