"""Real vascular integration, SI configuration and exact coefficient checkpoints."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .audit import sha256


def validate_reference(config):
    if config.get('condition_type')!='REFERENCE_NUMERICAL_CONDITION' or config.get('experimental') is not False:
        raise ValueError('Only the explicitly identified numerical reference condition is executable here')
    if config.get('warning')!='NOT_EXPERIMENTAL_PUMP_FLOW':
        raise ValueError('Missing numerical-reference flow provenance')
    physics=config['physics']
    vals=np.array([physics[k] for k in ('density_kg_m3','kinematic_viscosity_m2_s','dynamic_viscosity_pa_s','inlet_volume_flow_m3_s')],float)
    if not np.isfinite(vals).all() or np.any(vals<=0):raise ValueError('Positive finite SI properties and flow required')
    if not np.isclose(vals[2],vals[0]*vals[1],rtol=1e-14,atol=0):raise ValueError('Dynamic viscosity must equal rho * nu')
    bc=config['boundary_model']
    if bc['pressure_pin'] or bc['pressure_dirichlet'] or bc['inlet_velocity_profile'] is not None or bc['outlet_flow_split'] is not None:
        raise ValueError('Boundary model conflicts with the validated global-flow Stokes formulation')
    if len(bc['outlets'])!=3 or set(bc['outlets'].values())!={'sigma n=0'}:raise ValueError('Three natural zero-traction outlets required')
    return config


def load_config(root):
    root=Path(root);path=root/'configs/stage03_reference_vascular.yaml'
    snapshot=json.loads((root/'inputs/stage03/reference_condition.json').read_text())
    assert sha256(path)==snapshot['yaml_sha256'], 'Reference YAML differs from the frozen execution snapshot'
    config=validate_reference(snapshot['config'])
    assert sha256(root/'inputs/stage03/reference_source_config.yaml')==config['source']['configuration_sha256']
    return config,sha256(path)


def reynolds(config):
    p=config['physics'];g=config['inlet_geometry']
    mean=p['inlet_volume_flow_m3_s']/g['projected_area_m2']
    dh=4*g['projected_area_m2']/g['projected_perimeter_m']
    re=p['density_kg_m3']*mean*dh/p['dynamic_viscosity_pa_s']
    return {'inlet_area_m2':g['projected_area_m2'],'inlet_perimeter_m':g['projected_perimeter_m'],
            'hydraulic_diameter_m':dh,'mean_inlet_velocity_m_s':mean,'Re':re,
            'status':'STOKES_PHYSICS_REVIEW_REQUIRED' if re>=1 else 'LOW_RE_SUPPORTS_STOKES',
            'definition':'rho * (Q / projected inlet area) * (4 area / projected perimeter) / mu'}


def load_mesh(root,config,comm):
    from dolfinx import io,mesh
    from mpi4py import MPI
    root=Path(root);base=root/config['mesh']['source']
    for name,key in [('volume_mesh.npz','volume_mesh_sha256'),('fluid.xdmf','xdmf_sha256'),('fluid.h5','hdf5_sha256')]:
        assert sha256(base/'mesh'/name)==config['mesh'][key]
    with io.XDMFFile(comm,str(base/'mesh/fluid.xdmf'),'r') as f:
        domain=f.read_mesh(name='fluid_mesh',ghost_mode=mesh.GhostMode.shared_facet)
        domain.topology.create_connectivity(2,3)
        domain.topology.create_connectivity(3,2)
        ct=f.read_meshtags(domain,name='cell_tags');ft=f.read_meshtags(domain,name='facet_tags')
    nc=domain.topology.index_map(3).size_local;nf=domain.topology.index_map(2).size_local
    assert np.array_equal(ct.indices[ct.indices<nc],np.arange(nc)) and np.all(ct.values==100)
    exterior=mesh.exterior_facet_indices(domain.topology);exterior=exterior[exterior<nf]
    assert np.array_equal(np.sort(exterior),np.sort(ft.indices[ft.indices<nf]))
    counts={str(tag):comm.allreduce(int(np.count_nonzero((ft.values==tag)&(ft.indices<nf))),op=MPI.SUM) for tag in (1,2,3,4,5)}
    assert all(count>0 for count in counts.values()) and set(np.unique(ft.values))<={1,2,3,4,5}
    metadata={'gdim':domain.geometry.dim,'tdim':domain.topology.dim,'cells':domain.topology.index_map(3).size_global,
              'cell_tag':100,'boundary_counts':counts,'unlabelled_exterior_facets':0,'inlets':1,'outlets':3,
              'geometry_modified':False,'source_sha256':config['mesh']['volume_mesh_sha256']}
    return domain,ft,metadata


def save_primary(base,fields):
    """Gather owned nodal coefficients; no projection of the primary fields."""
    base=Path(base);comm=fields['velocity'].function_space.mesh.comm
    checkpoint={}
    for name in ('velocity','pressure'):
        f=fields[name];space=f.function_space
        n=space.dofmap.index_map.size_local;bs=space.dofmap.index_map_bs
        local=(space.tabulate_dof_coordinates()[:n].copy(),f.x.array[:n*bs].reshape(n,bs).copy())
        rows=comm.gather(local,root=0)
        if comm.rank==0:
            checkpoint[name+'_coordinates_m']=np.concatenate([r[0] for r in rows])
            checkpoint[name+'_values']=np.concatenate([r[1] for r in rows])
    path=base/'checkpoints/primary.npz'
    if comm.rank==0:
        assert not path.exists(),'Never overwrite a completed PDE checkpoint'
        path.parent.mkdir(parents=True,exist_ok=True)
        checkpoint['lambda_pa']=np.array(fields['lambda_pa'])
        np.savez_compressed(path,**checkpoint)
    comm.barrier()
    return sha256(path)


def restore_primary(root,config,comm):
    from dolfinx import fem
    from mpi4py import MPI
    from .spaces import create_spaces,real_owned_value
    root=Path(root);base=root/'outputs/stage03/reference'
    saved_meta=json.loads((base/'metadata/solve.json').read_text())
    path=base/'checkpoints/primary.npz'
    assert sha256(path)==saved_meta['checkpoint_sha256']
    assert sha256(root/'configs/stage03_reference_vascular.yaml')==saved_meta['config_sha256']
    domain,tags,mesh_meta=load_mesh(root,config,comm)
    spaces=create_spaces(domain);saved=np.load(path);fields={};distance_max={}
    for name,space in zip(('velocity','pressure'),spaces[:2]):
        n=space.dofmap.index_map.size_local;bs=space.dofmap.index_map_bs
        coordinates=space.tabulate_dof_coordinates()[:n]
        source=saved[name+'_coordinates_m']
        assert len(source)==space.dofmap.index_map.size_global and len(np.unique(source,axis=0))==len(source)
        d,indices=cKDTree(source).query(coordinates)
        tolerance=512*np.finfo(float).eps*np.max(np.abs(source))
        assert d.max(initial=0)<=tolerance and len(np.unique(indices))==len(indices)
        f=fem.Function(space,name=name);f.x.array[:n*bs]=saved[name+'_values'][indices].ravel();f.x.scatter_forward()
        assert np.array_equal(f.x.array[:n*bs],saved[name+'_values'][indices].ravel())
        fields[name]=f;distance_max[name]=comm.allreduce(float(d.max(initial=0)),op=MPI.MAX)
    multiplier=fem.Function(spaces[2],name='inlet_normal_traction_multiplier')
    multiplier.x.array[:spaces[2].dofmap.index_map.size_local]=float(saved['lambda_pa']);multiplier.x.scatter_forward()
    fields.update(multiplier=multiplier,lambda_pa=real_owned_value(multiplier))
    return domain,tags,fields,mesh_meta,distance_max
