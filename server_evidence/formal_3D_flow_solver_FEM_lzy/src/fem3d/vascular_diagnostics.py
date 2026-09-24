"""Solution integrals and physical-cell diagnostics for the real vascular mesh."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .diagnostics import integrator
from .derived_fields import derive
from .vascular import reynolds

PORT_TAGS={'inlet':4,'outlet_01':3,'outlet_02':5,'outlet_03':2}


def distribution(values,value):
    values=np.asarray(values,float)
    median=np.median(values,axis=0)
    ratio=np.divide(value,median,out=np.zeros_like(np.asarray(value,float)),where=median!=0)
    def plain(x):return np.asarray(x).tolist()
    return {'min':plain(np.min(values,axis=0)),'median':plain(median),'max':plain(np.max(values,axis=0)),
            'residual_over_neighbor_median':plain(ratio) if np.all(median!=0) else None,
            'ratio_defined':bool(np.all(median!=0))}


def local_monitor(root,domain,fields,derived,centers):
    """Find current owned cells from physical centroids, then use shared-face neighbors."""
    from mpi4py import MPI
    root=Path(root);comm=domain.comm;nc=domain.topology.index_map(3).size_local
    evidence=json.loads((root/'inputs/stage03/residual_source_evidence.json').read_text())
    targets=np.array([r['centroid_m'] for r in evidence['records']])
    distances,indices=cKDTree(centers[:nc]).query(targets)
    pairs=comm.allgather((distances,indices))
    owners=np.argmin(np.stack([v[0] for v in pairs]),axis=0)
    domain.topology.create_connectivity(3,2);domain.topology.create_connectivity(2,3)
    c2f=domain.topology.connectivity(3,2);f2c=domain.topology.connectivity(2,3)
    rows=[]
    for j,owner in enumerate(owners):
        if comm.rank!=owner:continue
        cell=int(indices[j]);distance=float(distances[j]);assert distance<=1e-15
        vertices=domain.geometry.x[domain.geometry.dofmaps[0][cell]]
        bary=np.linalg.solve((vertices[1:]-vertices[0]).T,targets[j]-vertices[0])
        assert np.all(bary>0) and bary.sum()<1
        neighbors=sorted({int(c) for f in c2f.links(cell) for c in f2c.links(f) if c!=cell})
        assert neighbors
        ids=np.array([cell]+neighbors,dtype=np.int32);points=centers[ids]
        u=fields['velocity'].eval(points,ids).reshape(-1,3)
        p=fields['pressure'].eval(points,ids).reshape(-1)
        grad=derived['velocity_gradient_s_inv'].eval(points,ids).reshape(-1,3,3)
        strain=derived['strain_rate_tensor_s_inv'].eval(points,ids).reshape(-1,3,3)
        vort=derived['vorticity_s_inv'].eval(points,ids).reshape(-1,3)
        vals={'velocity_m_s':u,'speed_m_s':np.linalg.norm(u,axis=1),'pressure_pa':p,
              'gradient_frobenius_s_inv':np.linalg.norm(grad,axis=(1,2)),
              'strain_frobenius_s_inv':np.linalg.norm(strain,axis=(1,2)),
              'vorticity_magnitude_s_inv':np.linalg.norm(vort,axis=1)}
        assert all(np.isfinite(v).all() for v in vals.values())
        rows.append({'source_centroid_m':targets[j].tolist(),'located_centroid_m':centers[cell].tolist(),
                     'centroid_match_displacement_m':distance,'source_quality':evidence['records'][j]['min_sicn'],
                     'owner_rank':int(owner),'current_local_cell':cell,'physical_vertex_coordinates_m':vertices.tolist(),
                     'reference_barycentric_coordinates':[float(1-bary.sum())]+bary.tolist(),
                     'neighbor_count':len(neighbors),'neighbor_centroids_m':centers[neighbors].tolist(),
                     'values':{k:np.asarray(v[0]).tolist() for k,v in vals.items()},
                     'neighbors':{k:distribution(v[1:],v[0]) for k,v in vals.items()},
                     'finite':True,'source_array_index_used_for_location':False})
    rows=[r for part in comm.allgather(rows) for r in part]
    rows.sort(key=lambda r:tuple(r['source_centroid_m']))
    assert len(rows)==len(targets)
    return {'status':'PASS','cells':rows,'count':len(rows),'all_finite':True,
            'location_method':'Current owned-cell physical centroid bijection plus positive reference barycentric membership; source array indices ignored',
            'neighborhood':'One shared-face ring, including valid MPI ghost neighbors',
            'manual_review':'MANUAL_RESIDUAL_CELL_REVIEW','automatic_ratio_threshold':None,
            'interpretation':'Finite local anomaly monitor only; one solution cannot prove absence of discretization error'}


def plane_pressure_sections(vertices,pressures,z_values,comm):
    """Exact P1 pressure mean on actual plane/tetra intersections (no mesh changes)."""
    from mpi4py import MPI
    edges=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)];rows=[]
    for z in z_values:
        crossing=np.flatnonzero((vertices[:,:,2].min(axis=1)<z)&(vertices[:,:,2].max(axis=1)>z))
        area_total=0.;integral=0.
        for cell in crossing:
            x=vertices[cell];p=pressures[cell];points=[];values=[]
            for a,b in edges:
                da=x[a,2]-z;db=x[b,2]-z
                if da*db<0:
                    t=-da/(db-da);points.append(x[a]+t*(x[b]-x[a]));values.append(p[a]+t*(p[b]-p[a]))
                elif da==0:
                    points.append(x[a]);values.append(p[a])
                elif db==0:
                    points.append(x[b]);values.append(p[b])
            points,ids=np.unique(np.array(points),axis=0,return_index=True);values=np.asarray(values)[ids]
            if len(points)<3:continue
            center=points.mean(axis=0);pc=values.mean();offset=points-center
            order=np.argsort(np.arctan2(offset[:,1],offset[:,0]));points=points[order];values=values[order]
            areas=np.linalg.norm(np.cross(points-center,np.roll(points,-1,axis=0)-center),axis=1)*.5
            area_total+=areas.sum();integral+=float(np.sum(areas*(values+np.roll(values,-1)+pc)/3))
        area_total=comm.allreduce(float(area_total),op=MPI.SUM);integral=comm.allreduce(integral,op=MPI.SUM)
        assert area_total>0 and np.isfinite([area_total,integral]).all()
        rows.append({'name':f'internal_z_{z*1e6:.3f}_um','z_m':float(z),'normal':[0,0,1],
                     'area_m2':area_total,'pressure_integral_pa_m2':integral,'mean_pressure_pa':integral/area_total,
                     'method':'Area integral of piecewise P1 pressure on all vascular intersections with this physical z plane',
                     'intersected_owned_cells':comm.allreduce(len(crossing),op=MPI.SUM)})
    return rows


def diagnose(root,config,domain,tags,fields):
    import ufl
    from dolfinx import fem
    from mpi4py import MPI
    root=Path(root);cache=root/'outputs/stage03/jit_cache';comm=domain.comm
    u,p=fields['velocity'],fields['pressure'];Q=config['physics']['inlet_volume_flow_m3_s']
    integrate=integrator(domain,cache)
    ds=ufl.Measure('ds',domain=domain,subdomain_data=tags,metadata={'quadrature_degree':6})
    dx=ufl.Measure('dx',domain=domain,metadata={'quadrature_degree':6});n=ufl.FacetNormal(domain)
    signed={name:integrate(ufl.dot(u,n)*ds(tag)) for name,tag in PORT_TAGS.items()}
    qin=-signed['inlet'];qout=sum(signed[name] for name in signed if name!='inlet')
    error=abs(qin-Q)/Q;closure=abs(qout-qin)/Q
    wall_flux=integrate(ufl.dot(u,n)*ds(1))
    fractions={name:signed[name]/qout for name in signed if name!='inlet'}
    flux={'Q_target_m3_s':Q,'Q_in_signed_m3_s':signed['inlet'],'Q_in_m3_s':qin,
          'Q_outlets_signed_m3_s':{k:v for k,v in signed.items() if k!='inlet'},'Q_out_total_m3_s':qout,
          'wall_flux_m3_s':wall_flux,'flow_fractions':fractions,
          'relative_inlet_error':error,'relative_mass_closure':closure,
          'status':'PASS' if max(error,closure)<=1e-10 else 'FAIL',
          'outlet_sign_review':'MANUAL_PHYSICS_REVIEW' if any(v<0 for k,v in signed.items() if k!='inlet') else 'ALL_NET_OUTFLOWS_POSITIVE',
          'method':'Computed P2 field integrated on physical tagged facets, quadrature degree 6 and MPI SUM'}
    V=u.function_space;wall_dofs=fem.locate_dofs_topological(V,2,tags.find(1));wall_dofs=wall_dofs[wall_dofs<V.dofmap.index_map.size_local]
    max_wall=comm.allreduce(float(np.abs(u.x.array.reshape(-1,3)[wall_dofs]).max(initial=0)),op=MPI.MAX)
    velocity_scale=Q/config['inlet_geometry']['projected_area_m2']
    finite=comm.allreduce(all(np.isfinite(f.x.array).all() for f in (u,p,fields['multiplier'])),op=MPI.LAND)
    grad_l2=float(np.sqrt(integrate(ufl.inner(ufl.grad(u),ufl.grad(u))*dx)))
    div_l2=float(np.sqrt(integrate(ufl.div(u)**2*dx)))
    boundary_pressure={}
    for name,tag in PORT_TAGS.items():
        area=integrate(1*ds(tag));value=integrate(p*ds(tag))
        boundary_pressure[name]={'area_m2':area,'pressure_integral_pa_m2':value,'mean_pressure_pa':value/area}
    derived,definitions=derive(u,cache)
    im=domain.topology.index_map(3);nc=im.size_local;nt=nc+im.num_ghosts
    vertices=domain.geometry.x[domain.geometry.dofmaps[0][:nt]]
    centers=vertices.mean(axis=1)
    cell_ids=np.arange(nc,dtype=np.int32)
    samples={'cell_centers_m':centers[:nc],'cell_velocity_m_s':u.eval(centers[:nc],cell_ids).reshape(-1,3),
             'cell_pressure_pa':p.eval(centers[:nc],cell_ids).reshape(-1)}
    for name,f in derived.items():samples[name]=f.eval(centers[:nc],cell_ids).reshape(nc,-1)
    samples_finite=comm.allreduce(all(np.isfinite(v).all() for v in samples.values()),op=MPI.LAND)
    monitor=local_monitor(root,domain,fields,derived,centers)
    bounds=np.array([domain.geometry.x[:,2].min(),domain.geometry.x[:,2].max()])
    zmin=comm.allreduce(float(bounds[0]),op=MPI.MIN);zmax=comm.allreduce(float(bounds[1]),op=MPI.MAX)
    zvalues=zmin+np.array([.27,.50,.73])*(zmax-zmin)
    vertex_p=p.eval(vertices[:nc].reshape(-1,3),np.repeat(cell_ids,4)).reshape(nc,4)
    sections=plane_pressure_sections(vertices[:nc],vertex_p,zvalues,comm)
    norms={'velocity_L2_norm_m_pow_2p5_s':float(np.sqrt(integrate(ufl.inner(u,u)*dx))),
           'pressure_integral_pa_m3':integrate(p*dx),'divergence_L2_m_pow_1p5_s':div_l2,
           'gradient_L2_m_pow_1p5_s':grad_l2,'relative_divergence_diagnostic':div_l2/(grad_l2+np.finfo(float).tiny)}
    scalars=[fields['lambda_pa'],max_wall,*norms.values(),*signed.values()]
    all_finite=bool(finite and samples_finite and np.isfinite(scalars).all())
    result={'status':'PASS' if all_finite and flux['status']=='PASS' and max_wall<=128*np.finfo(float).eps*velocity_scale else 'FAIL',
            'flux':flux,'lambda_pa':fields['lambda_pa'],'all_fields_finite':all_finite,
            'max_abs_wall_velocity_m_s':max_wall,'wall_roundoff_tolerance_m_s':128*np.finfo(float).eps*velocity_scale,
            'norms':norms,'boundary_pressure':boundary_pressure,'internal_pressure_sections':sections,
            'reynolds':reynolds(config),'residual_cells':monitor,'derived_definitions':definitions,
            'divergence_hard_gate':False,'wss_computed':False,'experimental':False}
    return result,samples


def export_visualization(base,fields,samples):
    """Export actual quadratic tetra connectivity and original P2 values for WSL rendering."""
    import basix.ufl
    from dolfinx import fem,plot
    base=Path(base);u=fields['velocity'];V=u.function_space;domain=V.mesh;comm=domain.comm
    scalar=fem.functionspace(domain,basix.ufl.element('Lagrange',domain.basix_cell(),2))
    p2=fem.Function(scalar);p2.interpolate(fields['pressure']);p2.x.scatter_forward()
    topology,types,points=plot.vtk_mesh(V)
    n=V.dofmap.index_map.size_local;ng=n+V.dofmap.index_map.num_ghosts;nc=domain.topology.index_map(3).size_local
    global_nodes=V.dofmap.index_map.local_to_global(np.arange(ng,dtype=np.int32))
    topology=topology.reshape(-1,11)[:nc];assert np.all(topology[:,0]==10)
    d,match=cKDTree(scalar.tabulate_dof_coordinates()).query(points[:n]);assert d.max(initial=0)<=1e-18
    local={'node_ids':global_nodes[:n],'points_m':points[:n],'velocity_m_s':u.x.array.reshape(-1,3)[:n],
           'pressure_pa':p2.x.array[match],'cells':global_nodes[topology[:,1:]],'cell_types':types[:nc],**samples}
    gathered=comm.gather(local,root=0)
    if comm.rank==0:
        count=V.dofmap.index_map.size_global
        output={name:np.empty((count,)+(gathered[0][name].shape[1:]),dtype=gathered[0][name].dtype)
                for name in ('points_m','velocity_m_s','pressure_pa')}
        for row in gathered:
            for name in output:output[name][row['node_ids']]=row[name]
        output['cells']=np.concatenate([row['cells'] for row in gathered])
        output['cell_types']=np.concatenate([row['cell_types'] for row in gathered])
        for name in samples:output[name]=np.concatenate([row[name] for row in gathered])
        assert all(np.isfinite(v).all() for v in output.values())
        (base/'solution').mkdir(parents=True,exist_ok=True)
        np.savez_compressed(base/'solution/fields_for_visualization.npz',**output)
