"""Independent evidence reader: never imports the converter or LAMMPS.
Checks actual population, actual data, dumps, execution logs and GPU activity.
All thresholds are read from the pre-run frozen stage contract.
"""
from pathlib import Path
import argparse,csv,hashlib,json,math,re,sqlite3
import numpy as np

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_json(p):return json.loads(Path(p).read_text())
def dump(p):
    with Path(p).open() as f:
        assert f.readline().strip()=='ITEM: TIMESTEP';step=int(f.readline())
        assert f.readline().strip()=='ITEM: NUMBER OF ATOMS';n=int(f.readline())
        assert f.readline().startswith('ITEM: BOX BOUNDS');box=np.array([[float(x) for x in f.readline().split()] for _ in range(3)])
        header=f.readline().split();assert header[:2]==['ITEM:','ATOMS'];names=header[2:]
        a=np.atleast_2d(np.loadtxt(f))
    assert len(a)==n and a.shape[1]==len(names)
    a=a[np.argsort(a[:,names.index('id')])]
    return {'step':step,'N':n,'names':names,'array':a,'box':box,'col':{name:a[:,i] for i,name in enumerate(names)}}
def vector(d,names):return np.column_stack([d['col'][k] for k in names])
def xyz(d):return vector(d,['xu','yu','zu'])
def velocity(d):return vector(d,['vx','vy','vz'])
def force(d):return vector(d,['fx','fy','fz'])
def maximum(a):return float(np.max(np.abs(a)))
def compare_dumps(a,b):
    return {'ids_exact':bool(np.array_equal(a['col']['id'],b['col']['id'])),'diameters_exact':bool(np.array_equal(a['col']['diameter'],b['col']['diameter'])),'max_position_difference_m':maximum(xyz(a)-xyz(b)),'max_velocity_difference_m_s':maximum(velocity(a)-velocity(b))}
def integrity(root):
    rows=[]
    for line in (root/'FROZEN_INPUT_SHA256SUMS').read_text().splitlines():
        expected,name=line.split('  ',1);p=root/name
        rows.append({'file':name,'pass':p.is_file() and sha(p)==expected})
    return {'status':'PASS' if all(r['pass'] for r in rows) else 'FAIL','files':rows}

def read_source(root,sp):
    with (root/sp['population']).open() as f:
        rows=list(csv.DictReader(f))
    ids=np.array([int(r['bubble_id'])+1 for r in rows]);d=np.array([float(r['diameter_um']) for r in rows])
    return ids,d

def initial_data(path):
    s=Path(path).read_text();a=s.split('Atoms # sphere\n\n')[1].split('\nVelocities')[0];v=s.split('\nVelocities\n\n')[1]
    atoms=np.atleast_2d(np.loadtxt(a.splitlines()));vel=np.atleast_2d(np.loadtxt(v.splitlines()))
    return atoms[np.argsort(atoms[:,0])],vel[np.argsort(vel[:,0])]

def kernel_summary(db):
    with sqlite3.connect('file:'+str(db)+'?mode=ro',uri=True) as c:
        tables={r[0] for r in c.execute("select name from sqlite_master where type='table'")}
        if 'CUPTI_ACTIVITY_KIND_KERNEL' not in tables:return {'status':'FAIL','kernel_count':0,'reason':'No CUDA kernel table'}
        cols={r[1] for r in c.execute('pragma table_info(CUPTI_ACTIVITY_KIND_KERNEL)')}
        name='demangledName' if 'demangledName' in cols else 'shortName'
        rows=c.execute(f'SELECT s.value, COUNT(*), SUM(k.end-k.start), MIN(k.start), MAX(k.end) FROM CUPTI_ACTIVITY_KIND_KERNEL k JOIN StringIds s ON s.id=k.{name} GROUP BY s.value ORDER BY SUM(k.end-k.start) DESC').fetchall()
    items=[{'name':r[0],'launches':r[1],'total_device_seconds':r[2]*1e-9,'first_ns':r[3],'last_ns':r[4]} for r in rows]
    nve=[r for r in items if 'NVESphere' in r['name']]
    gran=[r for r in items if 'GranHookeHistory' in r['name']]
    return {'status':'PASS' if nve and gran else 'FAIL','kernel_count':sum(r['launches'] for r in items),'nve_sphere_kernel_count':sum(r['launches'] for r in nve),'gran_hooke_history_kernel_count':sum(r['launches'] for r in gran),'kernels':items,'evidence_scope':'Case F complete launch including initialization and sparse outputs; no production performance claim'}

def evaluate_case(root,key):
    root=Path(root);ct=read_json(root/'contracts/MINIMAL_LAMMPS_PARTICLE_ENGINE_CONTRACT.json');g=ct['gates'];sp=next(x for x in ct['cases'] if x['case']==key);p=root/'cases'/key
    result={'case':key,'status':'FAIL','checks':{},'metrics':{}};checks=result['checks'];m=result['metrics']
    def check(k,v):checks[k]=bool(v)
    try:
        audit=integrity(root);check('frozen_inputs_unchanged',audit['status']=='PASS')
        a=dump(p/'initial.dump');b=dump(p/'final.dump');ids,dum=read_source(root,sp);data,veldata=initial_data(p/'particles.data')
        meta=read_json((root/sp['population']).with_suffix('.metadata.json'))
        popid=read_json(root/'provenance/POPULATION_SOURCE_IDENTITY.json')['populations'][sp['population_key']]
        sc=read_json(root/'provenance/frozen_sampler/contracts/SONOVUE_SAMPLER_CONTRACT_V0.json')
        check('population_provenance',sha(root/sp['population'])==meta['population_sha256']==popid['csv_sha256'] and meta['N']==sp['N'] and meta['seed']==sp['seed'] and meta['histogram_sha256']==sc['histogram_sha256'] and popid['sampler_contract_sha256']==ct['sampler_contract_sha256'])
        check('count_and_ids',a['N']==b['N']==len(ids)==sp['N'] and np.array_equal(a['col']['id'],ids) and np.array_equal(b['col']['id'],ids) and len(np.unique(ids))==len(ids))
        m['nonfinite_count']=int(np.sum(~np.isfinite(a['array']))+np.sum(~np.isfinite(b['array'])));check('nonfinite',m['nonfinite_count']==0)
        check('timestep',a['step']==0 and b['step']==sp['steps'])
        m['max_diameter_roundtrip_error_um']=maximum(a['col']['diameter']/1e-6-dum)
        check('diameter_roundtrip',m['max_diameter_roundtrip_error_um']<=g['diameter_roundtrip_um'])
        check('diameter_exactly_unchanged',np.array_equal(a['col']['diameter'],b['col']['diameter']))
        check('source_support',np.all((b['col']['diameter']/1e-6>=sc['support_min_um'])&(b['col']['diameter']/1e-6<=sc['support_max_um'])))
        check('radius_representation',np.array_equal(b['col']['radius']*2,b['col']['diameter']))
        expected_mass=sp['technical_density_kg_m3']*math.pi/6*b['col']['diameter']**3
        m['mass_max_relative_error']=maximum((b['col']['mass']-expected_mass)/expected_mass);check('mass_density',m['mass_max_relative_error']<=g['mass_relative_error'])
        check('actual_data_matches_source',np.array_equal(data[:,0],ids) and np.array_equal(data[:,2],dum*1e-6) and np.all(data[:,3]==sp['technical_density_kg_m3']))
        check('initial_dump_matches_actual_data',np.array_equal(xyz(a),data[:,4:7]) and np.array_equal(velocity(a),veldata[:,1:4]))
        check('omega',np.array_equal(vector(a,['omegax','omegay','omegaz']),veldata[:,4:7]) and np.array_equal(vector(a,['omegax','omegay','omegaz']),vector(b,['omegax','omegay','omegaz'])))
        inp=(p/'in.lammps').read_text();check('input_units_and_atom_style','units si\n' in inp and 'atom_style sphere\n' in inp)
        runtime=read_json(p/'RUN_METRICS.json');log=(p/'log.lammps').read_text()
        check('execution_complete',runtime['returncode']==0 and 'ENGINE_CASE_COMPLETED' in log and 'ERROR' not in log and not re.search(r'\b(?:nan|inf)\b',log,re.I))
        loop=re.findall(r'Loop time of ([\deE.+-]+) on (\d+) procs for (\d+) steps with (\d+) atoms',log)
        executed=[x for x in loop if int(x[2])==sp['steps'] and int(x[3])==sp['N']]
        check('actual_loop_steps_count',bool(executed))
        m['lost_atoms']=int(sp['N']-b['N']);check('lost_atoms',m['lost_atoms']==0 and 'Lost atoms:' not in log)
        m['neighbor_builds']=[int(x) for x in re.findall(r'Neighbor list builds = (\d+)',log)]
        m['dangerous_builds']=[int(x) for x in re.findall(r'Dangerous builds = (\d+)',log)]
        check('dangerous_builds',all(x==0 for x in m['dangerous_builds']))
        m['wall_seconds']=runtime['wall_seconds'];m['process_tree_peak_rss_kib']=runtime['process_tree_peak_rss_kib']
        if executed:
            m['solver_loop_seconds']=float(executed[-1][0]);m['solver_steps_per_s']=sp['steps']/m['solver_loop_seconds'] if m['solver_loop_seconds'] else None
        check('memory_recorded',m['process_tree_peak_rss_kib']>0)
        if sp['placement']=='grid':
            # Verify deterministic grid from actual data, including periodic spacing.
            shape=np.array(sp['grid_shape']);ijk=np.array(np.unravel_index(np.arange(len(ids)),tuple(shape))).T
            expected=(ijk+.5)*sp['spacing_m']
            check('actual_grid_positions',np.array_equal(data[:,4:7],expected))
            check('periodic_box',np.allclose(a['box'][:,1]-a['box'][:,0],np.array(shape)*sp['spacing_m'],rtol=0,atol=1e-18))
            m['minimum_periodic_center_spacing_m']=sp['spacing_m'];m['minimum_surface_gap_lower_bound_m']=sp['spacing_m']-float(np.max(dum))*1e-6
            check('initial_no_overlap',m['minimum_surface_gap_lower_bound_m']>sp['grid_clear_margin_m'])
            m['initial_overlap_count']=0 if checks['actual_grid_positions'] and checks['initial_no_overlap'] else None
        if not sp['contact']:
            target=data[:,4:7]+veldata[:,1:4]*sp['steps']*sp['dt_s']
            m['max_ballistic_position_error_m']=maximum(xyz(b)-target);m['max_velocity_error_m_s']=maximum(velocity(b)-veldata[:,1:4])
            check('ballistic_position',m['max_ballistic_position_error_m']<=g['ballistic_position_abs_m']);check('ballistic_velocity',m['max_velocity_error_m_s']<=g['ballistic_velocity_abs_m_s'])
            check('zero_force',maximum(force(a))<=g['stationary_force_abs_N'] and maximum(force(b))<=g['stationary_force_abs_N'])
        elif key!='case_c_technical_contact':
            check('stationary_position',maximum(xyz(b)-xyz(a))<=g['stationary_position_abs_m']);check('stationary_velocity',maximum(velocity(b))<=g['stationary_velocity_abs_m_s']);check('no_contact_force',maximum(force(a))<=g['stationary_force_abs_N'] and maximum(force(b))<=g['stationary_force_abs_N'])
        if sp['placement']=='two':
            rsum=.5*float(np.sum(a['col']['diameter']));d0=float(np.linalg.norm(xyz(a)[1]-xyz(a)[0]));d1=float(np.linalg.norm(xyz(b)[1]-xyz(b)[0]));m.update(initial_distance_m=d0,final_distance_m=d1,initial_overlap_m=max(0,rsum-d0),final_overlap_m=max(0,rsum-d1))
            if key=='case_b_two_separated':check('comfortable_gap',d0-rsum>9e-6)
            else:
                cc=read_json(root/'contracts/TECHNICAL_CONTACT_TEST_CONTRACT.json');f=force(a);fn=np.linalg.norm(f,axis=1);direction=xyz(a)[1]-xyz(a)[0]
                m['initial_forces_N']=f.tolist();m['force_symmetry_relative_error']=float(np.linalg.norm(f.sum(axis=0))/max(fn.max(),np.finfo(float).tiny))
                m['hooke_force_relative_error']=float(abs(fn[0]-cc['Kn_N_per_m']*m['initial_overlap_m'])/(cc['Kn_N_per_m']*m['initial_overlap_m']))
                check('controlled_overlap',m['initial_overlap_m']>0 and abs(m['initial_overlap_m']-cc['initial_overlap_m'])<1e-18)
                check('different_diameters',dum[0]!=dum[1]);check('finite_nonzero_initial_force',np.all(np.isfinite(f)) and np.all(fn>0))
                check('repulsive_direction',np.dot(f[0],direction)<0 and np.dot(f[1],direction)>0)
                check('force_symmetry',m['force_symmetry_relative_error']<=g['contact_force_symmetry_relative']);check('hooke_initial_force',m['hooke_force_relative_error']<=g['contact_initial_hooke_force_relative'])
                check('overlap_decreased',m['final_overlap_m']<m['initial_overlap_m'])
                vmax=cc['stability_preflight']['predicted_max_relative_speed_m_s']
                m['max_final_speed_m_s']=float(np.max(np.linalg.norm(velocity(b),axis=1)))
                check('no_explosion',m['max_final_speed_m_s']<=vmax*g['contact_speed_bound_multiple_of_energy_max'] and maximum(xyz(b)-xyz(a))<=g['contact_displacement_max_m'])
        if key.endswith('cpu_mpi4'):
            ref=dump(root/'cases/case_e_mpi_migration/cpu_mpi1/final.dump');cmp=compare_dumps(ref,b);m['mpi1_mpi4']=cmp
            check('mpi_agreement',cmp['ids_exact'] and cmp['diameters_exact'] and cmp['max_position_difference_m']<=g['mpi_position_abs_m'] and cmp['max_velocity_difference_m_s']<=g['mpi_velocity_abs_m_s'])
            m['ownership_changes']=int(np.sum(a['col']['proc']!=b['col']['proc']));m['initial_owners']={str(int(k)):int(v) for k,v in zip(*np.unique(a['col']['proc'],return_counts=True))};m['final_owners']={str(int(k)):int(v) for k,v in zip(*np.unique(b['col']['proc'],return_counts=True))}
            check('real_migration',len(m['initial_owners'])==4 and m['ownership_changes']>0 and '4 by 1 by 1 MPI processor grid' in log)
        if key=='cpu_gpu_comparison/gpu':
            ref=dump(root/'cases/cpu_gpu_comparison/cpu/final.dump');cmp=compare_dumps(ref,b);m['cpu_gpu']=cmp
            check('cpu_gpu_agreement',cmp['ids_exact'] and cmp['diameters_exact'] and cmp['max_position_difference_m']<=g['cpu_gpu_position_abs_m'] and cmp['max_velocity_difference_m_s']<=g['cpu_gpu_velocity_abs_m_s'])
        if sp['engine']=='gpu':check('kokkos_activated',re.search(r'KOKKOS mode(?: with Kokkos version [0-9.]+)? is enabled',log) is not None and re.search(r'(?:1|one) GPU',log,re.I) is not None)
        if key=='case_f_kokkos_gpu':
            ks=kernel_summary(p/'gpu_trace.sqlite');m['gpu_kernel_activity']={k:v for k,v in ks.items() if k!='kernels'};check('gpu_kernels_observed',ks['status']=='PASS')
            with (p/'RESOURCE_TRACE.csv').open() as f:tele=list(csv.DictReader(f))
            gpu=[r for r in tele if r['gpu_memory_mib'] and float(r['gpu_memory_mib'])>100]
            allgpu=[r for r in tele if r['gpu_memory_mib']]
            m['gpu_vram_peak_mib']=max(float(r['gpu_memory_mib']) for r in allgpu);m['gpu_utilization_peak_percent']=max(float(r['gpu_utilization_percent']) for r in allgpu)
            # Endpoints can include context setup/teardown. Use central half of GPU-resident samples.
            steady=gpu[len(gpu)//4:len(gpu)*3//4];mem=[float(r['gpu_memory_mib']) for r in steady];util=[float(r['gpu_utilization_percent']) for r in gpu]
            m['gpu_utilization_median_percent']=float(np.median(util)) if util else None;m['gpu_steady_samples']=len(mem)
            m['gpu_vram_steady_range_mib']=max(mem)-min(mem) if mem else None
            check('gpu_vram_stable',len(mem)>=g['gpu_vram_min_samples'] and m['gpu_vram_steady_range_mib']<=g['gpu_vram_post_initialization_range_mib_max'] and m['gpu_vram_peak_mib']<24564*g['gpu_vram_capacity_fraction_max'])
        result['status']='PASS' if all(checks.values()) else 'FAIL'
    except Exception as e:
        import traceback
        result['exception']=repr(e);result['traceback']=traceback.format_exc()
    return result

def evaluate_all(root):
    root=Path(root);ct=read_json(root/'contracts/MINIMAL_LAMMPS_PARTICLE_ENGINE_CONTRACT.json')
    cases=[evaluate_case(root,k) for k in ct['sequence']]
    build=read_json(root/'build_provenance/LAMMPS_BUILD_PROVENANCE.json')
    return {'status':'PASS' if all(c['status']=='PASS' for c in cases) and build['status']=='PASS' and integrity(root)['status']=='PASS' else 'FAIL','scope':'Independent numerical and recorded provenance audit; local transport SHA256 is a separate gate','cases':cases,'frozen_input_integrity':integrity(root)['status'],'build_provenance':build['status'],'nonfinite_count':sum(c['metrics'].get('nonfinite_count',0) for c in cases),'lost_atoms':sum(c['metrics'].get('lost_atoms',0) for c in cases)}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--case');ap.add_argument('--output',required=True);args=ap.parse_args()
    r=evaluate_case(args.root,args.case) if args.case else evaluate_all(args.root)
    Path(args.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps({'status':r['status'],'output':args.output}));raise SystemExit(0 if r['status']=='PASS' else 1)
