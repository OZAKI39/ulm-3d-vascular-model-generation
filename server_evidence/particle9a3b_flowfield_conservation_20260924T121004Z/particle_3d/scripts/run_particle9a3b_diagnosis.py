"""Six CPU slab audits, exact element divergence, native/export/Particle comparison."""
from pathlib import Path
import argparse
import csv
import json
import os
import platform
import time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

import numpy as np
import pyvista as pv

from particle_3d.audit import read_frozen
from particle_3d.field import FrozenFEMField
from particle_3d.interior_section import cut_tetrahedra, root_component, p1_flux
from particle_3d.flowfield_conservation_diagnosis import (p1_divergence,face_topology,boundary_owners,
    slab_balance,nearest_root_arclength,stats,clipped_boundary_flux,tetra_plane_polygon,polygon_flux)
from inventory_particle9a3b import write, checkpoint, sha

ENV = None


def csv_write(path, rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def section(grid,center,normal):
    sec,components,contains=root_component(cut_tetrahedra(grid,center,normal),center,1e-16)
    if not contains:raise ValueError('Root center absent from section')
    return p1_flux(sec,normal),components


def job(candidate):
    grid,field,div,vol,top,boundaries,representations,qin=ENV
    start=time.monotonic()
    center=np.array([float(candidate['center_m_'+i]) for i in 'xyz'])
    normal=np.array([float(candidate['normal_'+i]) for i in 'xyz'])
    s=float(candidate['arclength_m'])
    result,weights=slab_balance(field.points_m,field.tetra,field.velocity_nodes_m_s,div,vol,top,boundaries,center,normal)
    qwall=result['surface_fluxes']['WALL']['Q_m3_s'];qi=-result['surface_fluxes']['INLET']['Q_m3_s']
    unwanted=sum(abs(value['Q_m3_s']) for role,value in result['surface_fluxes'].items() if role.startswith('OUTLET'))
    unwanted_area=sum(value['area_m2'] for role,value in result['surface_fluxes'].items() if role.startswith('OUTLET'))
    rep=dict(section_s_m=s,candidate_id=int(candidate['candidate_id']),refinement=int(candidate['refinement']))
    for name,g in representations.items():
        value,components=section(g,center,normal)
        rep['Q_'+name+'_m3_s']=value['signed_Q_m3_s']
        rep['Q_positive_'+name+'_m3_s']=value['positive_Q_m3_s']
        rep['error_'+name]=abs(value['signed_Q_m3_s']/qin-1)
    qp=rep['Q_particle_m3_s'];D=result['volume_integral_divergence_m3_s']
    row=dict(section_s_m=s,candidate_id=rep['candidate_id'],refinement=rep['refinement'],
       Q_in_m3_s=qi,Q_section_m3_s=qp,Q_wall_m3_s=qwall,Q_loss_m3_s=qi-qp,
       volume_integral_divergence_m3_s=D,negative_volume_integral_divergence_m3_s=-D,
       loss_minus_wall_m3_s=qi-qp-qwall,
       closure_residual_m3_s=qp+qwall-qi-D,
       closure_residual_relative_inlet=(qp+qwall-qi-D)/qin,
       Q_loss_relative_inlet=(qi-qp)/qin,divergence_relative_inlet=D/qin,wall_relative_inlet=qwall/qin,
       Q_section_numpy_m3_s=result['Q_section_numpy_m3_s'],
       two_integrators_relative_difference=abs(qp-result['Q_section_numpy_m3_s'])/qin,
       inlet_reproduction_relative_difference=abs(qi-qin)/qin,other_outlet_area_m2=unwanted_area,
       other_outlet_abs_Q_m3_s=unwanted,clipped_volume_m3=result['clipped_volume_m3'],
       full_tetra_count=result['full_tetra_count'],partial_tetra_count=result['partial_tetra_count'],
       wall_triangle_count=result['surface_fluxes']['WALL']['triangle_count'],
       max_wall_normal_velocity_m_s=result['surface_fluxes']['WALL']['max_abs_normal_velocity_m_s'],
       worker_pid=os.getpid(),runtime_s=time.monotonic()-start)
    if unwanted_area!=0 or row['inlet_reproduction_relative_difference']>1e-10:
        raise ValueError('Slab closure contains an unexpected cap or incomplete inlet: '+str(row))
    if row['two_integrators_relative_difference']>1e-10 or abs(row['closure_residual_relative_inlet'])>1e-10:
        raise ValueError('Independent Gauss/plane integration failed: '+str(row))
    return row,rep,weights


def synthetic_actual_mesh(grid,field,top,boundaries,candidate,weights):
    """Reuse one of six real clipped volumes; three globally affine controls."""
    p,t=field.points_m,field.tetra
    c=np.array([float(candidate['center_m_'+i]) for i in 'xyz']);n=np.array([float(candidate['normal_'+i]) for i in 'xyz'])
    z=(p-p.mean(axis=0))*10. # [1/s] times position, thus m/s
    fields={'constant':(np.tile([.001,-.002,.0005],(len(p),1)),0.),
            'linear_div_free':(np.column_stack([z[:,0],-z[:,1],np.zeros(len(p))])+[.001,0,0],0.),
            'linear_nonzero_div':(z+[.001,0,0],30.)}
    centroids=p[t].mean(axis=1);mask=weights>0;dist=(p-c)@n
    crossing=np.flatnonzero(mask & (np.max(dist[t],axis=1)>0))
    rows=[]
    for name,(u,expected_div) in fields.items():
        div,vol=p1_divergence(p,t,u)
        synthetic=grid.copy();synthetic['Velocity']=u
        vtkq,_=section(synthetic,c,n)
        numpyq=sum(polygon_flux(*tetra_plane_polygon(p[t[i]],u[t[i]],c,n),n)[0] for i in crossing)
        boundary={role:clipped_boundary_flux(p,u,tri,owners,centroids,mask,c,n)['Q_m3_s']
                  for role,(tri,owners) in boundaries.items()}
        total=numpyq+sum(boundary.values());expected=expected_div*weights.sum()
        scale=max(abs(numpyq),sum(abs(q) for q in boundary.values()),abs(expected),1e-30)
        rows.append(dict(field=name,expected_divergence_s_inv=expected_div,max_element_divergence_error=float(abs(div-expected_div).max()),
          plane_flux_vtk=vtkq['signed_Q_m3_s'],plane_flux_numpy=numpyq,volume_divergence=float(weights@div),
          expected_volume_divergence=float(expected),boundary_outward_flux=total,
          analytic_gauss_residual_relative=float(abs(total-expected)/scale),
          integrators_relative_difference=float(abs(numpyq-vtkq['signed_Q_m3_s'])/scale),
          pass_check=bool(abs(total-expected)/scale<1e-10 and abs(numpyq-vtkq['signed_Q_m3_s'])/scale<1e-10)))
    if not all(x['pass_check'] for x in rows):raise ValueError('Synthetic actual-mesh check failed')
    return rows


def main(a):
    global ENV
    start=time.monotonic();report=a.repo/'particle_3d/reports/particle9a3b_flowfield_conservation';data=report/'data'
    assert (data/'velocity_representation_inventory.json').is_file(), 'Inventory must precede calculation'
    summary,grid,flow,surfaces=read_frozen(a.repo/'formal_3D_flow_solver/FEM_SimVascular')
    field=FrozenFEMField.from_grids(grid,flow);grid.points=field.points_m.copy();grid['Velocity']=field.velocity_nodes_m_s.copy()
    print('Inventory present. Computing exact divergence and tetra face connectivity.',flush=True)
    div,vol=p1_divergence(field.points_m,field.tetra,field.velocity_nodes_m_s)
    trace=np.trace(field.gradients_s_inv,axis1=1,axis2=2)
    write(data/'particle_gradient_crosscheck.json',dict(max_abs_divergence_difference=float(abs(div-trace).max()),
        method='independent batched solve of edge equations vs FrozenFEMField inverse/einsum',finite_difference=False))
    top=face_topology(field.tetra);boundaries={}
    for role,surface in surfaces.items():
        ids=np.asarray(surface['GlobalNodeID'],int)-1;tri=ids[surface.faces.reshape(-1,4)[:,1:]]
        boundaries[role]=(tri,boundary_owners(top,tri))
    tri,owners=boundaries['INLET'];xyz=field.points_m[tri];av=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])/2
    centroids=field.points_m[field.tetra].mean(axis=1)
    flip=np.einsum('ij,ij->i',av,xyz.mean(axis=1)-centroids[owners])<0;av[flip]*=-1
    uv=field.velocity_nodes_m_s[tri];qin=float(-np.einsum('ij,ij->',av,uv.mean(axis=1)))
    inlet_area=float(np.linalg.norm(av,axis=1).sum());inward=-av.sum(axis=0)/np.linalg.norm(av.sum(axis=0))
    saved=json.loads((a.repo/'formal_3D_flow_solver/FEM_SimVascular/frozen_reference/new_flow_flux.json').read_text())['integrated_inlet_Q_m3_s']
    cp,h=checkpoint(a.case/'run/1-procs/stFile_071.bin',len(field.points_m))
    native=pv.read(a.case/'run/1-procs/result_071.vtu');native.points=np.asarray(native.points,float)
    exported=pv.read(a.case/'frozen_flow/steady_flow_mean_2p0_mmps.vtu');exported.points=np.asarray(exported.points,float)
    nativecp=native.copy();nativecp['Velocity']=cp[:,:3].copy()
    representations={'native_checkpoint':nativecp,'native_vtu':native,'export':exported,'particle':grid}
    write(data/'inlet_flux_audit.json',dict(inlet_area_m2=inlet_area,Q_in_m3_s=qin,saved_Q_in_m3_s=saved,
          relative_saved_difference=abs(qin-saved)/saved,mean_normal_velocity_m_s=qin/inlet_area,
          max_speed_m_s=float(np.linalg.norm(uv,axis=2).max()),inward_normal=inward.tolist(),
          inlet_surface_point_arrays=list(surfaces['INLET'].point_data.keys()),
          point_data_velocity_note='INLET mesh cap has identifiers, no independent Velocity array; values mapped by GlobalNodeID from same volume nodal DOFs',
          nodal_velocity_all_representations_equal=all(np.array_equal(g['Velocity'],field.velocity_nodes_m_s) for g in representations.values())))
    walltri,wallowner=boundaries['WALL'];wx=field.points_m[walltri];wu=field.velocity_nodes_m_s[walltri]
    wav=np.cross(wx[:,1]-wx[:,0],wx[:,2]-wx[:,0])/2
    wav[np.einsum('ij,ij->i',wav,wx.mean(axis=1)-centroids[wallowner])<0]*=-1
    warea=np.linalg.norm(wav,axis=1);wn=wav/warea[:,None];wq=np.einsum('ij,ij->i',wav,wu.mean(axis=1))
    wallnodes=np.unique(walltri);wspeed=np.linalg.norm(field.velocity_nodes_m_s[wallnodes],axis=1)
    csv_write(data/'wall_node_velocity.csv',[dict(global_node_id=int(i)+1,
        velocity_x_m_s=field.velocity_nodes_m_s[i,0],velocity_y_m_s=field.velocity_nodes_m_s[i,1],
        velocity_z_m_s=field.velocity_nodes_m_s[i,2]) for i in wallnodes])
    write(data/'wall_flux_audit.json',dict(node_count=len(wallnodes),triangle_count=len(walltri),
          all_wall_velocity_components_exactly_zero=bool(np.all(wu==0)),max_wall_node_speed_m_s=float(wspeed.max()),
          max_wall_normal_velocity_m_s=float(abs(np.einsum('tvi,ti->tv',wu,wn)).max()),
          total_wall_flux_m3_s=float(wq.sum()),sum_absolute_triangle_flux_m3_s=float(abs(wq).sum()),
          max_absolute_triangle_flux_m3_s=float(abs(wq).max()),normal_orientation='outward from adjacent tetra centroid'))
    csv_write(data/'wall_flux_audit.csv',[dict(wall_triangle_id=i,tetra_id=int(wallowner[i]),area_m2=warea[i],Q_wall_m3_s=wq[i],
          max_abs_normal_velocity_m_s=float(abs(wu[i]@wn[i]).max())) for i in range(len(walltri))])
    walladj=np.zeros(len(div),bool);walladj[wallowner]=True
    root=json.loads((a.p9a3/'data/root_topology.json').read_text());rootpoints=np.array(root['points_m'])
    arc,rootdist=nearest_root_arclength(centroids,rootpoints)
    daughterdist=np.full(len(div),np.inf)
    for x,y in np.array(root['daughter_segments_m']):
        e=y-x;t=np.clip((centroids-x)@e/np.dot(e,e),0,1)
        daughterdist=np.minimum(daughterdist,np.linalg.norm(centroids-x-t[:,None]*e,axis=1))
    rootmask=(rootdist<=daughterdist) & (arc<root['arclength_m'][-1])
    records=dict(domain=stats(div,vol),root_inlet_to_first_junction=stats(div[rootmask],vol[rootmask]),
        wall_adjacent=stats(div[walladj],vol[walladj]),interior=stats(div[~walladj],vol[~walladj]),
        root_wall_adjacent=stats(div[rootmask&walladj],vol[rootmask&walladj]),
        root_interior=stats(div[rootmask&~walladj],vol[rootmask&~walladj]),
        root_assignment='tetra centroid nearest authoritative root polyline rather than any daughter edge; arclength before junction. Descriptive anatomical bins, NOT the volume mask used for Gauss.',
        wall_adjacent_definition='tetra with at least one actual WALL triangle; not merely a wall node')
    rows=list(csv.DictReader((a.p9a3/'data/interior_section_search.csv').open()))
    legal=[r for r in rows if r['geometry_pass']=='True']
    targets=[.07e-6,1e-6,10e-6,30e-6,50e-6,max(float(r['arclength_m']) for r in legal)]
    chosen=[min(legal,key=lambda r:abs(float(r['arclength_m'])-x)) for x in targets]
    write(data/'selected_candidates.json',dict(rule='nearest geometrically legal old candidate to 0.07,1,10,30,50 um; last legal before junction; selected without inspecting new flux',candidates=chosen))
    ENV=(grid,field,div,vol,top,boundaries,representations,qin)
    with ProcessPoolExecutor(max_workers=a.workers,mp_context=multiprocessing.get_context('fork')) as pool:
        results=list(pool.map(job,chosen))
    balance=[x[0] for x in results];reps=[x[1] for x in results]
    csv_write(data/'slab_gauss_balance.csv',balance);csv_write(data/'internal_flux_by_representation.csv',reps)
    lastweights=results[-1][2];pick=lastweights>0
    records['exact_inlet_to_last_legal_section_slab']=stats(div[pick],lastweights[pick])
    write(data/'divergence_statistics.json',records)
    # Streaming CSV avoids hundreds of thousands of dictionaries in memory.
    with (data/'tetra_divergence.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['tetra_id','centroid_x_m','centroid_y_m','centroid_z_m','volume_m3','divergence_s_inv',
           'root_arclength_m','distance_to_root_m','wall_adjacent','root_segment','last_legal_slab_volume_fraction'])
        w.writerows((i,*centroids[i],vol[i],div[i],arc[i],rootdist[i],int(walladj[i]),int(rootmask[i]),lastweights[i]/vol[i]) for i in range(len(div)))
    bins=[]
    for lo,hi in zip(np.linspace(0,root['arclength_m'][-1],57)[:-1],np.linspace(0,root['arclength_m'][-1],57)[1:]):
        keep=rootmask&(arc>=lo)&(arc<hi)
        if keep.any():
            quant=np.quantile(div[keep],[.01,.1,.5,.9,.99]);bins.append(dict(s_um=(lo+hi)*.5e6,count=int(keep.sum()),
                P01=quant[0],P10=quant[1],median=quant[2],P90=quant[3],P99=quant[4],volume_weighted_mean=float(vol[keep]@div[keep]/vol[keep].sum())))
    csv_write(data/'root_divergence_bins.csv',bins)
    print('All six real slabs close. Running affine controls on the same real mesh.',flush=True)
    controls=synthetic_actual_mesh(grid,field,top,boundaries,chosen[3],results[3][2])
    write(data/'synthetic_real_mesh_validation.json',controls)
    # Whole-domain identity is an independent check of cell orientation and boundary roles.
    allq={}
    for role,(tri,owners) in boundaries.items():
        x=field.points_m[tri];av=np.cross(x[:,1]-x[:,0],x[:,2]-x[:,0])/2
        av[np.einsum('ij,ij->i',av,x.mean(axis=1)-centroids[owners])<0]*=-1
        allq[role]=float(np.einsum('ij,ij->',av,field.velocity_nodes_m_s[tri].mean(axis=1)))
    write(data/'global_gauss_balance.json',dict(signed_outward_boundary_fluxes=allq,
        volume_integral_divergence_m3_s=float(vol@div),closure_residual_m3_s=float(sum(allq.values())-vol@div),
        global_inlet_outlet_relative_error=abs(sum(allq.values()))/qin))
    manifest=dict(status='NATIVE_FIELD_LOCAL_CONSERVATION_LIMITATION_FOUND',Q_in_m3_s=qin,
        max_gauss_residual_relative=max(abs(r['closure_residual_relative_inlet']) for r in balance),
        max_two_integrators_difference=max(r['two_integrators_relative_difference'] for r in balance),
        representation_field_identity=True,wall_exactly_zero=bool(np.all(wu==0)),
        section_flux_loss_relative=[r['Q_loss_relative_inlet'] for r in balance],workers=a.workers,
        actual_worker_pids=sorted({r['worker_pid'] for r in balance}),thread_environment={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']},
        runtime_s=time.monotonic()-start,host=platform.node(),python=platform.python_version(),numpy=np.__version__,pyvista=pv.__version__,
        source_repo=str(a.repo.resolve()),case_path=str(a.case.resolve()),utc_unix_s=time.time(),
        tolerance_recommendation='NO_TOLERANCE_RECOMMENDATION_YET',injection_status='BLOCKED_BY_NATIVE_FIELD_CONSERVATION',
        production_changes=False,bubble_events_generated=0,bubble_trajectories_run=0)
    write(data/'diagnosis_summary.json',manifest);print(json.dumps(manifest,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--case',type=Path,required=True)
    p.add_argument('--p9a3',type=Path,required=True);p.add_argument('--workers',type=int,default=6);main(p.parse_args())
