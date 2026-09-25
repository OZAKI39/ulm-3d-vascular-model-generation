"""One-time WSL preparation, before geometry selection and before numerical cases."""
from pathlib import Path
import numpy as np,h5py,vtk,json,hashlib,sys,shutil,datetime,subprocess
from vtk.util.numpy_support import vtk_to_numpy
P=Path(__file__).resolve().parents[1];SC=Path('/home/lzy/projects/sonovue_size_distribution_v0');SM=Path('/home/lzy/projects/compre_output/pure_fluid_new_medium_smoke/20260915_161937');RV=Path('/home/lzy/projects/compre_output/review_sync/20260913-013144/publish_checkout/review_bundle/step3_review')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,x):p.parent.mkdir(exist_ok=True,parents=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
for root in [SC,Path('/home/lzy/projects/compre_output/lammps_particle_engine/20260915_214759')]:
 r=subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=root,text=True,capture_output=True);assert r.returncode==0,r.stderr
sys.path.insert(0,str(SC/'src'));from sonovue_sampler import SonoVueDistribution,write_population
sc=json.loads((SC/'contracts/SONOVUE_SAMPLER_CONTRACT_V0.json').read_text())
# Freeze identity/selection before examining geometry fit; no redraw or resizing.
selection={'A':{'seed':20260918,'N':1,'rule':'first sample'},'B':{'seed':20260922,'N':1,'rule':'first sample'},'C':{'seed':20260919,'N':100,'rule':'minimum absolute diameter distance to frozen target d50; tie => smallest bubble_id'},'D':{'seed':20260920,'N':10000,'rule':'minimum absolute distance to frozen target d10, d50, d90; tie => smallest bubble_id'},'E':{'seed':20260921,'N':1000,'rule':'all samples'}}
put(P/'contracts/PARTICLE_SELECTION_CONTRACT.json',{'status':'FROZEN_BEFORE_GEOMETRY_FIT','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'rules':selection,'histogram_sha256':sc['histogram_sha256'],'sampler_contract_sha256':sha(SC/'contracts/SONOVUE_SAMPLER_CONTRACT_V0.json')})
chosen={}
for name,info in selection.items():
 p=P/'populations'/(name+'.csv');write_population(SonoVueDistribution(),info['N'],info['seed'],p,role='SAMPLER_VALIDATION_ONLY')
 a=np.genfromtxt(p,names=True,delimiter=',');a=np.atleast_1d(a)
 if name in ['A','B']:idx=[0]
 elif name=='C':idx=[int(np.argmin(abs(a['diameter_um']-sc['target_quantiles_um']['d50'])))]
 elif name=='D':idx=[int(np.argmin(abs(a['diameter_um']-sc['target_quantiles_um'][q]))) for q in ['d10','d50','d90']]
 else:idx=list(range(len(a)))
 chosen[name]={'seed':info['seed'],'source_N':info['N'],'source_bubble_ids':idx,'diameters_um':a['diameter_um'][idx].tolist(),'population_sha256':sha(p),'sampler_contract_sha256':sha(SC/'contracts/SONOVUE_SAMPLER_CONTRACT_V0.json'),'histogram_sha256':sc['histogram_sha256'],'role':'COUPLING_ENGINEERING_TEST_ONLY'}
put(P/'provenance/POPULATION_IDENTITY.json',chosen)
for name in ['contracts/SONOVUE_SAMPLER_CONTRACT_V0.json','src/sonovue_sampler.py','input/FROZEN_SONOVUE_HISTOGRAM.csv']:
 out=P/'provenance/frozen_sampler'/name;out.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(SC/name,out)
num=json.loads((SM/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json').read_text());geom=json.loads((RV/'step2/diagnostics/palabos_geometry.json').read_text());origin=np.array(geom['physical_origin_m']);dx=geom['effective_dx_m'];dims=np.array(geom['lattice_shape'],dtype=np.int64)
for source,name in [(SM/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json','PBS_BSA_NUMERICS.json'),(RV/'step2/diagnostics/palabos_geometry.json','PALABOS_GEOMETRY.json'),(RV/'step2/physical_to_lattice_transform.md','PHYSICAL_TO_LATTICE_TRANSFORM.md'),(SM/'source/mpi_support.hpp','ORIGINAL_SNAPSHOT_WRITER.hpp'),(SM/'source/gpu_adapter.hpp','ORIGINAL_GPU_ADAPTER.hpp'),(SM/'source/vascularPoC.cpp','ORIGINAL_PALABOS_DRIVER.cpp')]:shutil.copyfile(source,P/'provenance'/name)
stl=RV/'step1/geometry_contract/geometry/cfd_surface_axis_aligned_inlet_m.stl';shutil.copyfile(stl,P/'provenance/closed_geometry_m.stl')
mask_source=RV/'step2/closed_flag_matrix.vti';reader=vtk.vtkXMLImageDataReader();reader.SetFileName(str(mask_source));reader.Update();im=reader.GetOutput();assert np.array_equal(np.array(im.GetDimensions())-1,dims)
assert im.GetCellData().GetNumberOfArrays()==1
mask=vtk_to_numpy(im.GetCellData().GetArray(0)).reshape(-1);fluid=np.flatnonzero(mask).astype(np.uint64);assert len(fluid)==182694
source=SM/'run_5000/diagnostics/field_samples/fields_5000.bin'
with source.open('rb') as f:
 n=int(np.fromfile(f,dtype='<u8',count=1)[0]);raw=np.fromfile(f,dtype=[('id','<u8'),('rho','<f8'),('u','<f8',(3,))])
assert n==len(raw)==len(fluid) and np.array_equal(raw['id'],fluid)
assert np.all(np.isfinite(raw['u']))
u=raw['u']*(num['dx_m']/num['dt_s'])
# No use of the invalid derived magnitude; physical vector is recovered directly.
metadata={'format_version':'FROZEN_FLOW_FIELD_V0','source_case':'PBS_BSA_5000_STEP_SMOKE','source_file':str(source),'source_sha256':sha(source),'vector_source':'Velocity_m_s','vector_recovery':'PackedNode.u[3] in frozen binary snapshot, multiplied once by frozen dx_m/dt_s; no derived magnitude input','coordinate_units':'m','velocity_units':'m/s','dx_m':dx,'origin_m':origin.tolist(),'nx':int(dims[0]),'ny':int(dims[1]),'nz':int(dims[2]),'coordinate_convention':'global integer Palabos node: x_m = physical_origin_m + (ix,iy,iz)*dx_m; no half-cell shift','flattening':'linear_index = ix + nx*(iy + ny*iz); x fastest','iteration':5000,'physical_time_s':5000*num['dt_s'],'medium':'1x PBS + 1% BSA at 25 C, development assumption, not experimentally measured','scientific_status':'ENGINEERING_TRANSIENT_FIELD_ONLY','NEW_MEDIUM_INLET_MULTIPLIER_VALIDATION':'NOT_PERFORMED','TRANSIENT_BACKFLOW_OBSERVED':'YES','geometry_sha256':sha(stl),'mask_source':str(mask_source),'mask_source_sha256':sha(mask_source),'PBS_BSA_numerics_sha256':sha(SM/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json'),'Palabos_GPU_commit':'4127697e90169bbef982295f1d1c933cf6e90caa','missing_policy':'in-bounds index absent from fluid_linear_index = SOLID; present as fluid but absent from linear_index = MISSING','source_population_nodes':len(fluid),'representation':'sparse structured, complete fluid classification, no dense zero-velocity array'}
def write_field(path,indices,vel,fluid_ids,dims,origin,dx,meta):
 with h5py.File(path,'x') as h:
  h.attrs['format_version']='FROZEN_FLOW_FIELD_V0';h.attrs['coordinate_units']='m';h.attrs['velocity_units']='m/s';h.attrs['metadata_json']=json.dumps(meta)
  h.create_dataset('dims',data=np.asarray(dims,dtype='<i8'));h.create_dataset('origin_m',data=np.asarray(origin,dtype='<f8'));h.create_dataset('dx_m',data=float(dx))
  h.create_dataset('linear_index',data=np.asarray(indices,dtype='<u8'),compression='gzip',shuffle=True)
  h.create_dataset('fluid_linear_index',data=np.asarray(fluid_ids,dtype='<u8'),compression='gzip',shuffle=True)
  h.create_dataset('Velocity_m_s',data=np.asarray(vel,dtype='<f8'),compression='gzip',shuffle=True)
real=P/'fields/FROZEN_FLOW_FIELD_V0.h5';write_field(real,raw['id'],u,fluid,dims,origin,dx,metadata)
put(P/'contracts/FROZEN_FLOW_FIELD_CONTRACT.json',dict(metadata,field_file='fields/FROZEN_FLOW_FIELD_V0.h5',field_sha256=sha(real),status='FROZEN'))
rng=np.random.default_rng(2026091800);pick=rng.choice(fluid,size=256,replace=False);ijk=np.column_stack((pick%dims[0],(pick//dims[0])%dims[1],pick//(dims[0]*dims[1])));mapped=origin+ijk*dx;actual=[]
for idx in pick:
 bounds=im.GetCell(int(idx)).GetBounds();actual.append([(bounds[2*j]+bounds[2*j+1])/2 for j in range(3)])
err=float(np.max(np.abs(np.asarray(actual)-mapped)));assert err<=1e-12
np.savetxt(P/'validation/COORDINATE_NODE_CHECK.csv',np.column_stack((pick,ijk,mapped,np.asarray(actual))),delimiter=',',header='linear_index,ix,iy,iz,mapped_x,mapped_y,mapped_z,vtk_x,vtk_y,vtk_z',comments='',fmt=['%d']*4+['%.17g']*6)
put(P/'validation/COORDINATE_MAPPING_AUDIT.json',{'status':'PASS','reference':'Original Step2 CellData mask VTI cell centers; independently generated from frozen physical transform, shared nodes verified identical to original 5000-step snapshot','sampled_nodes':256,'seed':2026091800,'max_coordinate_error_m':err,'original_vti_origin_m':im.GetOrigin(),'original_vti_spacing_m':im.GetSpacing(),'lattice_node_origin_m':origin.tolist(),'field_node_indices_equal_Step2_closed_fluid_mask':True,'fluid_count':len(fluid),'half_cell_shift_applied_to_particle_mapping':False})
# Synthetic fields use the identical HDF5 schema and same C++ sampler path.
sdims=np.array([21,21,21]);sdx=8e-6;so=np.array([-80e-6]*3);si=np.arange(np.prod(sdims),dtype=np.uint64);xyz=so+np.column_stack((si%21,(si//21)%21,si//441))*sdx
intercept=np.array([1e-3,-2e-4,3e-4]);matrix=np.array([[20000.,10000.,-5000.],[-8000.,15000.,7000.],[3000.,-4000.,-12000.]])
synthetic={'uniform_velocity_m_s':intercept.tolist(),'linear_intercept_m_s':intercept.tolist(),'linear_matrix_per_s':matrix.tolist(),'dims':sdims.tolist(),'origin_m':so.tolist(),'dx_m':sdx,'status':'FROZEN_BEFORE_TESTS'}
for name,vel in [('UNIFORM_FIELD',np.tile(intercept,(len(si),1))),('LINEAR_FIELD',intercept+xyz@matrix.T)]:
 write_field(P/'fields'/(name+'.h5'),si,vel,si,sdims,so,sdx,dict(synthetic,source_case=name,scientific_status='SYNTHETIC_ENGINEERING_TEST_ONLY'))
missing=int(10+21*(10+21*10));keep=si!=missing
write_field(P/'fields/MISSING_CORNER_FIELD.h5',si[keep],np.tile(intercept,(np.sum(keep),1)),si,sdims,so,sdx,dict(synthetic,source_case='SYNTHETIC_MISSING_CORNER'))
put(P/'contracts/SYNTHETIC_FIELDS_CONTRACT.json',dict(synthetic,files={p.name:sha(p) for p in (P/'fields').glob('*FIELD.h5')},missing_corner_index=missing))
rho=1000.;mu=.001;tmin=rho*(sc['support_min_um']*1e-6)**2/(18*mu);tmax=rho*(sc['support_max_um']*1e-6)**2/(18*mu);dt=tmin/500
put(P/'contracts/PARTICLE_TIMESTEP_CONTRACT.json',{'status':'FROZEN_BEFORE_FORMAL_RUN','technical_density_kg_m3':rho,'mu_pa_s':mu,'support_um':[sc['support_min_um'],sc['support_max_um']],'tau_min_s':tmin,'tau_max_s':tmax,'dt_particle_s':dt,'definition':'tau_min/500, more conservative than tau_min/50; accommodates velocity-dependent force evaluation in standard LAMMPS nve/sphere without changing force law','dt_not_taken_from_LBM':True,'dt_LBM_s':num['dt_s'],'normalization_velocity':'norm(uniform fluid velocity) or norm(initial linear fluid velocity)','normalization_position':'characteristic_velocity * particle_tau, not absolute coordinate norm'})
print(json.dumps({'FIELD_EXPORT':'PASS','field_bytes':real.stat().st_size,'coordinate_error_m':err,'C_chosen_diameter_um':chosen['C']['diameters_um'][0],'tau_min_s':tmin,'tau_max_s':tmax,'particle_dt_s':dt},indent=2))
