"""Run once in the validated WSL NumPy environment, before any engine case."""
from pathlib import Path
import sys,json,hashlib,shutil,datetime,math,subprocess
import numpy as np
P=Path(__file__).resolve().parents[1]
S=Path('/home/lzy/projects/sonovue_size_distribution_v0')
sys.path.insert(0,str(S/'src'));sys.path.insert(0,str(P/'src'))
from sonovue_sampler import SonoVueDistribution,write_population
from prepare_lammps_particles import write_case,load_population

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,data):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(data,indent=2)+'\n')

assert sys.version_info[:2]==(3,12) and np.__version__=='1.26.4'
subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=S,check=True,capture_output=True)
sc=json.loads((S/'contracts/SONOVUE_SAMPLER_CONTRACT_V0.json').read_text())
for name in ['src/sonovue_sampler.py','contracts/SONOVUE_SAMPLER_CONTRACT_V0.json','input/FROZEN_SONOVUE_HISTOGRAM.csv','validation/SONOVUE_SAMPLER_VALIDATION.json']:
 out=P/'provenance/frozen_sampler'/name;out.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(S/name,out)
pops={}; dist=SonoVueDistribution()
for name,n,seed in [('single',1,20260915),('pair',2,20260916),('thousand',1000,20260915),('gpu',50000,20260917)]:
 out=P/'populations'/(name+'.csv');m=write_population(dist,n,seed,out,role='SAMPLER_VALIDATION_ONLY')
 pops[name]={'file':'populations/'+name+'.csv','N':n,'seed':seed,'csv_sha256':sha(out),'metadata_sha256':sha(out.with_suffix('.metadata.json')),'sampler_contract_sha256':sha(S/'contracts/SONOVUE_SAMPLER_CONTRACT_V0.json'),'histogram_sha256':sc['histogram_sha256'],'sampler_source_sha256':sha(S/'src/sonovue_sampler.py'),'engine_population_role':'GPU_ENGINE_SMOKE_ONLY' if n==50000 else 'ENGINE_TEST_ONLY'}
put(P/'provenance/POPULATION_SOURCE_IDENTITY.json',{'status':'PASS','sampling_executed_on':'original WSL environment','python':sys.version,'numpy':np.__version__,'sampler_API':'SonoVueDistribution + write_population','no_new_sampling_logic':True,'sampler_files_unmodified':True,'populations':pops})
_,pair=load_population(P/pops['pair']['file']);mass=1000*math.pi/6*pair**3;meff=mass.prod()/mass.sum();kn=.001;dt=1e-8;overlap=.02*pair.min();omega=math.sqrt(kn/meff)
contact={'contract_name':'TECHNICAL_CONTACT_TEST_CONTRACT','status':'FROZEN_BEFORE_RUN','pair_style':'gran/hooke/history','gpu_counterpart':'gran/hooke/history/kk','Kn_N_per_m':kn,'Kt_N_per_m':0,'gamma_n_per_s':0,'gamma_t_per_s':0,'friction_coefficient':0,'dampflag':0,'limit_damping':False,'initial_overlap_m':float(overlap),'initial_overlap_definition':'0.02 * minimum of the two sampled diameters','dt_s':dt,'steps':200,'source_document':'build_provenance/pinned_source_excerpts/doc/src/pair_gran.rst','interpretation':'Undamped, frictionless, normal Hooke repulsion. Technical machinery test only.','stability_preflight':{'mass_kg':mass.tolist(),'effective_mass_kg':float(meff),'omega_per_s':omega,'dt_times_omega':dt*omega,'gate_dt_times_omega_max':.05,'status':'PASS' if dt*omega<.05 else 'FAIL','verlet_linear_oscillator_stability_limit':2,'predicted_max_relative_speed_m_s':overlap*omega},'BUBBLE_BUBBLE_PHYSICS_VALIDATED':'NO','SONOVUE_CONTACT_MATERIAL_PROPERTIES':'NOT_DEFINED'}
assert contact['stability_preflight']['status']=='PASS'
put(P/'contracts/TECHNICAL_CONTACT_TEST_CONTRACT.json',contact)
put(P/'contracts/LAMMPS_PARTICLE_UNIT_CONTRACT.json',{'contract_name':'LAMMPS_PARTICLE_UNIT_CONTRACT','status':'FROZEN_BEFORE_RUN','units':'si','length':'m','time':'s','mass':'kg','force':'N','velocity':'m/s','density':'kg/m3','diameter_conversion':'diameter_um * 1e-6','atom_style':'sphere','mass_definition':'density * pi / 6 * diameter_m**3','TECHNICAL_TEST_DENSITY_KG_M3':1000,'TECHNICAL_TEST_DENSITY_IS_MICROBUBBLE_PHYSICS':'NO','fixed_radius':True,'dynamic_radius':False,'shape_deformation':False,'id_mapping':'particle_id = bubble_id + 1'})

def spec(key,pop,placement,steps,dt_s=1e-8,contact=True,ranks=1,vel=None,shape=None,gap=None):
 sp={'case':key,'population':pops[pop]['file'],'population_key':pop,'N':pops[pop]['N'],'seed':pops[pop]['seed'],'steps':steps,'dt_s':dt_s,'contact':contact,'mpi_ranks':ranks,'placement':placement,'technical_density_kg_m3':1000,'velocity_m_s':vel or [0.,0.,0.],'omega_rad_s':[0.,0.,0.],'box_m':[100e-6]*3,'periodic':[True]*3,'LAMMPS_UNITS':'si','ROLE':pops[pop]['engine_population_role']}
 if shape:sp.update(grid_shape=shape,spacing_m=8e-6,grid_clear_margin_m=1e-6,box_m=(np.array(shape)*8e-6).tolist())
 if gap is not None:sp['surface_gap_m']=gap
 return sp
specs=[spec('diameter_roundtrip','thousand','grid',0,shape=[10]*3),
 spec('case_a_single_ballistic','single','single',200,1e-6,False,vel=[.001,-.002,.003]),
 spec('case_b_two_separated','pair','two',200,gap=10e-6),
 spec('case_c_technical_contact','pair','two',200,gap=-float(overlap)),
 spec('case_d_polydisperse_1000','thousand','grid',200,shape=[10]*3),
 spec('case_e_mpi_migration/cpu_mpi1','thousand','grid',1000,1e-5,False,1,[.004,.001,-.002],[10]*3),
 spec('case_e_mpi_migration/cpu_mpi4','thousand','grid',1000,1e-5,False,4,[.004,.001,-.002],[10]*3),
 spec('case_f_kokkos_gpu','gpu','grid',1000,shape=[50,50,20]),
 spec('cpu_gpu_comparison/cpu','thousand','grid',200,1e-5,False,1,[.004,.001,-.002],[10]*3),
 spec('cpu_gpu_comparison/gpu','thousand','grid',200,1e-5,False,1,[.004,.001,-.002],[10]*3)]
for sp in specs:
 sp['engine']='gpu' if sp['case'] in ['case_f_kokkos_gpu','cpu_gpu_comparison/gpu'] else 'cpu'
 receipt=write_case(P/sp['population'],sp,P/'cases'/sp['case']);put(P/'cases'/sp['case']/'CONVERSION_RECEIPT.json',receipt)
gates={'diameter_roundtrip_um':1e-10,'ballistic_position_abs_m':1e-10,'ballistic_velocity_abs_m_s':1e-10,'stationary_position_abs_m':1e-14,'stationary_velocity_abs_m_s':1e-14,'stationary_force_abs_N':1e-20,'contact_force_symmetry_relative':1e-10,'contact_initial_hooke_force_relative':1e-10,'contact_speed_bound_multiple_of_energy_max':1.1,'contact_displacement_max_m':10e-6,'mass_relative_error':1e-12,'mpi_position_abs_m':1e-10,'mpi_velocity_abs_m_s':1e-10,'cpu_gpu_position_abs_m':1e-10,'cpu_gpu_velocity_abs_m_s':1e-9,'lost_atoms':0,'nonfinite':0,'gpu_vram_post_initialization_range_mib_max':256,'gpu_vram_min_samples':3,'gpu_vram_capacity_fraction_max':.9,'total_output_bytes_max':1024**3,'run_wall_time_limit_s':180,'nsys_export_time_limit_s':120}
metadata=json.loads((P.parent/'OFFICIAL_RELEASE_METADATA.json').read_text())
put(P/'contracts/MINIMAL_LAMMPS_PARTICLE_ENGINE_CONTRACT.json',{'contract_name':'MINIMAL_LAMMPS_PARTICLE_ENGINE_V0','status':'FROZEN_BEFORE_FIRST_FORMAL_RUN','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'LAMMPS_source':metadata,'unit_contract':'LAMMPS_PARTICLE_UNIT_CONTRACT.json','technical_contact_contract':'TECHNICAL_CONTACT_TEST_CONTRACT.json','sampler_contract_sha256':sha(S/'contracts/SONOVUE_SAMPLER_CONTRACT_V0.json'),'technical_density_kg_m3':1000,'sequence':[s['case'] for s in specs],'cases':specs,'gates':gates,'neighbor_policy':'every 1 delay 0 check no; skin 1e-6 m','noninteracting_pair':'zero cutoff 6e-6 m; force identically zero','gpu_flags':['-k','on','g','1','-sf','kk','-pk','kokkos','neigh','half','newton','off','gpu/aware','off'],'float_precision':'double','mpi_process_grid':'mpi_ranks x 1 x 1','resource_sampling_seconds':.1,'gpu_activity_evidence':'single Case F execution under existing Nsight Systems; CUDA kernel activity table','automatic_case_retries':False,'failed_foundation_stops_dependents':True,'cases_are_not_scientific_concentration':True,'Palabos':'NOT_USED','LATBOLTZ':'OFF','fluid_drag':'OFF','wall':'OFF','RBC':'OFF','adhesion':'OFF','performance_optimization':'NO'})
print(json.dumps({'PREPARATION':'PASS','contact_stability':contact['stability_preflight'],'population_count':len(pops),'case_count':len(specs)},indent=2))
