#!/usr/bin/env python3
"""Derive the reference condition from the current read-only 3D source config."""
import json
import shutil
import sys
from pathlib import Path
import numpy as np
import yaml
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json

R,I=ROOT/'reports/stage03',ROOT/'inputs/stage03'
assert (R/'pytest_after_cleanup.xml').is_file()
assert json.loads((R/'cleanup_audit.json').read_text())['status']=='PASS'
source=Path('/home/lzy/projects/ulm_3D_vascular/configs/cfd_flow.yaml')
current=yaml.safe_load(source.read_text())
lineage=json.loads((ROOT/'reports/stage00/source_contract.json').read_text())
source_run=(source.parent.parent/current['paths']['source_surface_run']).resolve()
assert source_run==Path(lineage['source_run']).resolve()
rho=float(current['physics']['density_kg_m3'])
nu=float(current['physics']['kinematic_viscosity_m2_s'])
Q=float(current['boundary_conditions']['target_volume_flow_m3_s'])
assert min(rho,nu,Q)>0 and np.isfinite([rho,nu,Q]).all()
assert np.isclose(rho*Q,float(current['boundary_conditions']['target_mass_flow_kg_s']),rtol=1e-12,atol=0)
mesh=ROOT/'outputs/stage01_7/selected'
ports=json.loads((mesh/'planar_port_contract_v2.json').read_text())
inlet=ports['ports']['inlet']
rim=np.array(inlet['rim_coordinates_m']);basis=np.array(inlet['basis'])
xy=(rim-np.array(inlet['plane_origin_m']))@basis[:2].T
area=abs(float(np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1))/2))
assert np.isclose(area,inlet['formal_projected_area_m2'],rtol=1e-12,atol=0)
perimeter=float(np.linalg.norm(np.roll(xy,-1,axis=0)-xy,axis=1).sum())
config={'schema_version':1,'name':'stage03_reference_vascular','condition_type':'REFERENCE_NUMERICAL_CONDITION',
        'experimental':False,'warning':'NOT_EXPERIMENTAL_PUMP_FLOW',
        'physics':{'density_kg_m3':rho,'kinematic_viscosity_m2_s':nu,'dynamic_viscosity_pa_s':rho*nu,
                   'inlet_volume_flow_m3_s':Q,'dynamic_viscosity_definition':'rho * nu'},
        'source':{'configuration_path':str(source),'configuration_sha256':sha256(source),
                  'density_key':'physics.density_kg_m3','kinematic_viscosity_key':'physics.kinematic_viscosity_m2_s',
                  'volume_flow_key':'boundary_conditions.target_volume_flow_m3_s',
                  'stage00_lineage_path':'reports/stage00/source_contract.json',
                  'stage00_lineage_sha256':sha256(ROOT/'reports/stage00/source_contract.json'),
                  'geometry_source_run':str(source_run),'outlet_pressures_imported':False},
        'mesh':{'source':'outputs/stage01_7/selected','volume_mesh_sha256':sha256(mesh/'mesh/volume_mesh.npz'),
                'xdmf_sha256':sha256(mesh/'mesh/fluid.xdmf'),'hdf5_sha256':sha256(mesh/'mesh/fluid.h5'),
                'port_contract_sha256':sha256(mesh/'planar_port_contract_v2.json')},
        'inlet_geometry':{'projected_area_m2':area,'projected_perimeter_m':perimeter,
                          'hydraulic_diameter_m':4*area/perimeter,'algebraic_length_scale_m':float(np.sqrt(area/np.pi)),
                          'length_scale_definition':'Equivalent inlet radius from projected area; algebraic scaling only'},
        'boundary_model':{'wall':'u=0','inlet':'integral(u dot outward n)=-Q; sigma n=-lambda n',
                          'outlets':{n:'sigma n=0' for n in ('outlet_01','outlet_02','outlet_03')},
                          'inlet_velocity_profile':None,'outlet_flow_split':None,'pressure_pin':False,
                          'pressure_dirichlet':False},
        'solver':{'mpi_ranks':4,'omp_num_threads':1,'ksp_type':'preonly','pc_type':'lu','factor_solver':'mumps',
                  'gpu_used':False,'unvalidated_iterative_fallback':False},
        'gates':{'relative_inlet_error_max':1e-10,'relative_mass_closure_max':1e-10,'local_finiteness_only':True,
                 'divergence_hard_threshold':None,'outlet_positive_flux_hard_gate':False},
        'created_utc':timestamp()}
path=ROOT/'configs/stage03_reference_vascular.yaml'
assert not path.exists()
path.write_text(yaml.safe_dump(config,sort_keys=False,allow_unicode=True))
write_json(I/'reference_condition.json',{'yaml_sha256':sha256(path),'config':config})
template={'schema_version':1,'name':'stage03_experiment_template','condition_type':'EXPERIMENTAL',
          'experimental':True,'inlet_volume_flow_m3_s':None,
          'instruction':'User must explicitly supply measured syringe-pump flow; no reference-flow fallback.',
          'status':'NOT_EXECUTED'}
(ROOT/'configs/stage03_experiment_template.yaml').write_text(yaml.safe_dump(template,sort_keys=False))
shutil.copy2(source,I/'reference_source_config.yaml')
write_json(R/'reference_condition_provenance.json',{'status':'PASS','timestamp':timestamp(),'source':config['source'],
    'physics':config['physics'],'inlet_geometry':config['inlet_geometry'],'config_sha256':sha256(path),
    'experimental_run':'NOT_EXECUTED','source_mass_flow_consistent':True})
q=json.loads((mesh/'qc/volume_quality.json').read_text())['quality']
residual=[r for r in q['worst_elements'] if r['min_sicn']<.1]
write_json(I/'residual_source_evidence.json',{'source':'outputs/stage01_7/selected/qc/volume_quality.json',
    'source_sha256':sha256(mesh/'qc/volume_quality.json'),'threshold':.1,'records':residual})
print(json.dumps({'physics':config['physics'],'inlet_geometry':config['inlet_geometry'],'residual_evidence':residual},indent=2))
