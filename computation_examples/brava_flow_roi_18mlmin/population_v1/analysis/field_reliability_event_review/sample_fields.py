"""Production GPU FEM query, independent strict containment and closed CV flux."""
from common import *
import gc
import itertools
import pyvista as pv
import vtk
import torch
from particle_3d.field import FrozenFEMField
from particle_3d.population_v1.gpu_field import GPUField
from particle_3d.population_v1.gpu_wall import GPUWall
from particle_3d.wall_geometry import WallGeometry
from particle_3d.convex_triangle import triangle_closest_many
from particle_3d.population_v1.bvh import tetra_adjacency
from prepare_coordinates import expand
from bounds_index import BoundsIndex

def extended_weights(field,x,ids):
    """Cancellation control on unmoved SI coordinates, only for strict ties."""
    t=np.asarray(field.points_m[field.tetra[ids]],np.longdouble)
    v=np.asarray(x,np.longdouble)-t[:,0];e=t[:,1:]-t[:,:1]
    det=np.sum(e[:,0]*np.cross(e[:,1],e[:,2]),axis=1)
    tail=np.stack([np.sum(v*np.cross(e[:,1],e[:,2]),axis=1),
                   np.sum(e[:,0]*np.cross(v,e[:,2]),axis=1),
                   np.sum(e[:,0]*np.cross(e[:,1],v),axis=1)],axis=1)/det[:,None]
    return np.c_[1-tail.sum(1),tail]

def strict_ids(field,x,ids):
    ids=ids.copy();valid=ids>=0
    ix=np.flatnonzero(valid);w=field.geometry.weights(x[ix],ids[ix])
    bad=ix[(w.min(1)<0)|(w.max(1)>1)]
    # A tolerance-first production choice can lie just outside one cell while
    # another cell strictly contains the original, unmoved query point.
    for i in bad:
        pad=field.geometry.candidate_padding_m; candidates=vtk.vtkIdList()
        field._locator.FindCellsWithinBounds(np.c_[x[i]-pad,x[i]+pad].ravel(),candidates)
        cs=np.array(sorted(candidates.GetId(j) for j in range(candidates.GetNumberOfIds())),int)
        ws=extended_weights(field,x[i],cs);ok=np.flatnonzero((ws.min(1)>=0)&(ws.max(1)<=1))
        ids[i]=int(cs[ok[0]]) if len(ok) else -1
    return ids

def sample(field,gpu,x):
    allids=[];allflow=[]
    for start in range(0,len(x),65536):
        xx=np.asarray(x[start:start+65536],float)
        result=gpu.sample(torch.tensor(xx,device='cuda',dtype=torch.float64))
        original=result[0].cpu().numpy();data=result[1].cpu().numpy()
        ids=strict_ids(field,xx,original)
        changed=np.flatnonzero((ids>=0)&(ids!=original))
        if len(changed):
            weights=field.geometry.weights(xx[changed],ids[changed]);nodes=field.tetra[ids[changed]]
            data[changed,:3]=np.einsum('ni,nij->nj',weights,field.velocity_nodes_m_s[nodes])
            data[changed,3]=np.einsum('ni,ni->n',weights,field.pressure_nodes_pa[nodes])
            data[changed,4:13]=field.gradients_s_inv[ids[changed]].reshape(-1,9)
            data[changed,13:16]=field.curl_s_inv[ids[changed]]
            data[changed,16:25]=field.strain_s_inv[ids[changed]].reshape(-1,9)
        data[ids<0]=np.nan
        allids.append(ids);allflow.append(data)
    return np.concatenate(allids),np.concatenate(allflow)

def wall_samples(wall,gpu,x,detail=False):
    tid,values,*_=gpu.query(torch.tensor(x,device='cuda',dtype=torch.float64),torch.zeros(len(x),device='cuda',dtype=torch.float64))
    tid=tid.cpu().numpy();values=values.cpu().numpy()
    feature=(values[:,10:13]>128*np.finfo(float).eps).sum(1)
    # Edge/vertex finite-feature normals are used by the existing geometry;
    # the distinct facet normals are recorded, never mistaken for flow reversal.
    ties=[];ambiguous=np.zeros(len(x),bool)
    if detail:
        for i,p in enumerate(x):
            d=values[i,0];pad=wall.roundoff_m;c=vtk.vtkIdList()
            wall._locator.FindCellsWithinBounds(np.c_[p-d-pad,p+d+pad].ravel(),c)
            ids=np.array(sorted(c.GetId(j) for j in range(c.GetNumberOfIds())),int)
            cp,bary=triangle_closest_many(p,wall.triangles[ids]);dist=norm(cp-p)
            near=ids[dist<=dist.min()+pad]
            different=np.max(norm(wall.normal_in[near]-wall.normal_in[near[0]]))>128*np.finfo(float).eps
            ambiguous[i]=different
            if len(near)>1:ties.append(dict(index=i,triangles=near.tolist(),inward_facet_normals=wall.normal_in[near].tolist(),different_facet_normals=bool(different)))
    else:ambiguous=feature<3
    return dict(triangle=tid,d_m=values[:,0],nearest_m=values[:,1:4],normal=values[:,7:10],barycentric=values[:,10:13],feature=feature,normal_ambiguous=ambiguous),ties

def make_wall():
    poly=pv.read(paths('M0')[0]/'SV_MESH/mesh-surfaces/WALL.vtp');faces=poly.faces.reshape(-1,4)[:,1:]
    wall=WallGeometry(poly.points[faces],provenance={'role':'original_BraVa_WALL'})
    wall.global_node_ids=poly.point_data['GlobalNodeID'][faces]-1
    return wall,GPUWall(wall)

def coordinates(pid,kind):
    if kind=='accepted': return np.load(OLD/f'original_{pid}.npz')['x_m']
    return np.load(HERE/f'{kind}_{pid}.npz')['x_m']

def normals_stage():
    wall,gpu=make_wall();protocol=read(HERE/'EVENT_FIELD_PROTOCOL.json');ties={}
    for pid in PIDS:
        for kind in ['accepted','path','plane','profiles','rise']:
            x=coordinates(pid,kind)
            values,t=wall_samples(wall,gpu,x,detail=kind in ['accepted','profiles'])
            if kind=='profiles':
                e=protocol['events'][str(pid)];profile=np.load(HERE/f'profiles_{pid}.npz');valid=np.ones(len(x),bool)
                for k in [0,1]:
                    ix=np.flatnonzero(profile['profile_index']==k);n=np.array(e['normal_profiles']['normals'][k]);o=np.array(e['normal_profiles']['origins_m'][k])
                    valid[ix]=(values['feature'][ix]==3)&(norm(values['nearest_m'][ix]-o)<=8*wall.roundoff_m)
                    valid[ix[0]]=True
                    valid[ix]=np.logical_and.accumulate(valid[ix])
                values['original_normal_ray_valid']=valid
            np.savez_compressed(HERE/f'wall_{kind}_{pid}.npz',**values)
            if t:ties[f'{pid}_{kind}']=t
            print('WALL',pid,kind,len(x),len(t),flush=True)
    put('WALL_NORMAL_TIES.json',dict(roundoff_m=wall.roundoff_m,features='3 FACE,2 EDGE,1 VERTEX; exact original geometry',
       tie_scope='All accepted states and independent normal profiles; interior quadrature/plane/rise retain feature ambiguity flag; tie list not exhaustively enumerated there.',ties=ties))

def triangle_rule(level):
    n=2**level;rows=[]
    for i in range(n):
        for j in range(n-i):
            rows.append([(i+1/3)/n,(j+1/3)/n])
            if i+j<n-1:rows.append([(i+2/3)/n,(j+2/3)/n])
    tail=np.array(rows);return np.c_[1-tail.sum(1),tail]

def fluxes(field,gpu,pid,layer):
    z=np.load(HERE/f'cv_{pid}_L{layer}.npz');tri=z['triangles_m'];normals=z['normal_out'];areas=z['area_m2'];wall=z['wall_face']
    levels=[]
    # M0 faces are exactly linear and their native vertex IDs are known.
    # All physical boundary facets are also identical across admitted meshes;
    # use their exact nodal boundary trace instead of evaluating rounded points
    # just outside a zero-thickness cap. No interior query is projected.
    exact_faces=np.full(len(tri),np.nan)
    if getattr(field,'review_name',None)=='M0':
        exact_faces=np.einsum('ij,ij->i',field.velocity_nodes_m_s[z['original_node_ids']].mean(1),normals)*areas
    else:
        from scipy.spatial import cKDTree
        tree=cKDTree(field.points_m)
        for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
            poly=pv.read(paths('M0')[0]/f'SV_MESH/mesh-surfaces/{role}.vtp')
            keys={tuple(sorted(row)) for row in (poly.point_data['GlobalNodeID'][poly.faces.reshape(-1,4)[:,1:]]-1)}
            for fi,nodeids in enumerate(z['original_node_ids']):
                if tuple(sorted(nodeids)) not in keys:continue
                distance,ids=tree.query(tri[fi]);assert np.all(distance==0)
                exact_faces[fi]=float(field.velocity_nodes_m_s[ids].mean(0)@normals[fi]*areas[fi])
    exact_faces[wall]=0
    for level in range(7):
        bary=triangle_rule(level);q=exact_faces.copy();invalid=0;numeric_floor=0.
        for fi in np.flatnonzero(~np.isfinite(exact_faces)):
            x=bary@tri[fi];ids,f=sample(field,gpu,x);bad=ids<0;invalid+=int(bad.sum())
            if np.any(bad):q[fi]=np.nan;continue
            integrand=f[:,:3]@normals[fi]
            q[fi]=areas[fi]*integrand.mean()
            numeric_floor+=areas[fi]*np.sum(np.abs(integrand))*np.finfo(float).eps*32
        levels.append(dict(level=level,points_per_face=len(bary),invalid=invalid,per_face_flux_m3_s=q.tolist(),net_flux_m3_s=float(q.sum()),
            throughput_m3_s=float(.5*np.abs(q).sum()),roundoff_scale_m3_s=float(numeric_floor)))
    return dict(event=pid,layer=layer,volume_m3=float(z['volume_m3']),tetra=len(z['cells']),faces=len(tri),wall_faces=int(wall.sum()),exact_trace_faces=int(np.isfinite(exact_faces).sum()),levels=levels,
      unavailable_reason=None if levels[-1]['invalid']==0 else 'Unresolved strict containment in surface quadrature; retained, not extrapolated',
      final_net_refinement_difference_m3_s=levels[-1]['net_flux_m3_s']-levels[-2]['net_flux_m3_s'],
      final_face_L1_refinement_difference_m3_s=float(np.abs(np.array(levels[-1]['per_face_flux_m3_s'])-levels[-2]['per_face_flux_m3_s']).sum()))

def field_stage(name):
    comp=read(HERE/'FIELD_COMPARABILITY.json');meta=comp['fields'][name]
    mesh=pv.read(meshpath(name));flow=pv.read(flowpath(name));native=pv.read(BASE/meta['raw_payload']['native'])
    # Same original GlobalNodeID correspondence. Never use nearest nodes.
    if 'GlobalNodeID' in native.point_data:
        gids=np.asarray(native.point_data['GlobalNodeID'],int);order=np.argsort(gids)
        assert np.array_equal(gids[order],np.arange(1,mesh.n_points+1))
    else:
        # Historical native/frozen output omits GlobalNodeID. The established
        # pair validator checks exact unique coordinates and tetra node sets.
        from particle_3d.audit import validate_pair
        validate_pair(mesh,native)
        order=np.arange(mesh.n_points)
    assert np.array_equal(native.points[order],mesh.points)
    assert np.array_equal(native.point_data['Velocity'][order],flow.point_data['Velocity'])
    assert np.array_equal(native.point_data['Pressure'][order],flow.point_data['Pressure'])
    del native
    field=FrozenFEMField.from_grids(mesh,flow);del flow
    field._locator=BoundsIndex(field.points_m,field.tetra)
    field.review_name=name
    gpu=GPUField(field);torch.cuda.reset_peak_memory_stats()
    mesh.point_data['u_m_s']=field.velocity_nodes_m_s
    mesh.cell_data['G_s_inv']=field.gradients_s_inv.reshape(-1,9)
    mesh.cell_data['div_s_inv']=np.trace(field.gradients_s_inv,axis1=1,axis2=2)
    mesh.cell_data['curl_s_inv']=field.curl_s_inv
    mesh.cell_data['canonical_cell']=np.arange(mesh.n_cells)
    edges=np.array(list(itertools.combinations(range(4),2)))
    lengths=norm(field.points_m[field.tetra[:,edges[:,1]]]-field.points_m[field.tetra[:,edges[:,0]]])
    altitudes=1/norm(field.geometry.shape_gradients_m_inv)
    metrics=np.c_[lengths.min(1),lengths.max(1),altitudes.min(1),altitudes.max(1),field.geometry.volume_m3]
    results={};balance=[];adj=tetra_adjacency(field.tetra)
    for pid in PIDS:
        cells=[];coverage={}
        for kind in ['accepted','path','plane','profiles','rise']:
            x=coordinates(pid,kind);ids,data=sample(field,gpu,x);valid=ids>=0
            cellmetrics=np.full((len(x),5),np.nan);cellmetrics[valid]=metrics[ids[valid]]
            np.savez_compressed(HERE/f'{name}_{kind}_{pid}.npz',cell=ids,valid=valid,flow=data,cell_geometry=cellmetrics)
            cells.extend(ids[valid]);coverage[kind]=dict(total=len(x),valid=int(valid.sum()),invalid=int((~valid).sum()),invalid_indices=np.flatnonzero(~valid).tolist(),reason='OUTSIDE_REFERENCE_DOMAIN or strict boundary roundoff: no extrapolation; exact wall trace handled separately')
            print('SAMPLED',name,pid,kind,len(x),int((~valid).sum()),flush=True)
        roi=expand(np.unique(cells),adj,2);mesh.extract_cells(roi).save(HERE/f'vtk/{name}_sampled_roi_{pid}.vtu')
        results[pid]=coverage
        for layer in [0,1,2]:
            b=fluxes(field,gpu,pid,layer)
            if name=='M0':
                cv=np.load(HERE/f'cv_{pid}_L{layer}.npz');nodes=cv['original_node_ids']
                exact=np.einsum('ij,ij->i',field.velocity_nodes_m_s[nodes].mean(1),cv['normal_out'])*cv['area_m2']
                volume=np.sum(np.trace(field.gradients_s_inv[cv['cells']],axis1=1,axis2=2)*field.geometry.volume_m3[cv['cells']])
                b['M0_exact_per_face_m3_s']=exact.tolist();b['M0_exact_net_m3_s']=float(exact.sum());b['M0_div_volume_m3_s']=float(volume)
                b['M0_identity_error_m3_s']=float(exact.sum()-volume)
            balance.append(b);print('CV',name,pid,layer,b['levels'][-1]['net_flux_m3_s'],b['final_face_L1_refinement_difference_m3_s'],flush=True)
    peak=torch.cuda.max_memory_allocated()/2**30
    assert peak<8
    put(f'{name}_SAMPLING.json',dict(native_node_values_exact=True,canonical_pair_validated=True,coverage=results,GPU_peak_allocated_GiB=peak,
       geometry_columns=['min_edge_m','max_edge_m','min_altitude_m','max_altitude_m','volume_m3'],
       field_columns='u(0:3),p(3),G row-major(4:13),curl(13:16),strain(16:25)',
       strict_outside_policy=True,units='SI',velocity_roundoff_scale_m_s=float(128*np.finfo(float).eps*np.max(norm(field.velocity_nodes_m_s)))))
    put(f'{name}_BALANCE.json',finite_json(balance))
    # This is the native P1 negative regression only, not a new acceleration target.
    if name=='M0':
        put('M0_NEGATIVE_REGRESSION.json',dict(cell=138584,all_four_nodal_u_zero=bool(np.all(field.velocity_nodes_m_s[field.tetra[138584]]==0)),
          cell388769_G=field.gradients_s_inv[388769].tolist(),cell388769_div=float(np.trace(field.gradients_s_inv[388769])),
          cell388769_volume_m3=float(field.geometry.volume_m3[388769])))
    print('FIELD COMPLETE',name,peak,flush=True)

if __name__=='__main__':
    if sys.argv[1]=='normals':normals_stage()
    else:field_stage(sys.argv[1])
