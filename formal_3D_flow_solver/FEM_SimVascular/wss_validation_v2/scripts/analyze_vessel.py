"""Raw production WSS, fixed physical regions, exact P1 sections and adjacency."""
import argparse,csv,json,subprocess,sys
from pathlib import Path
import numpy as np
import pyvista as pv
from case_common import *
from vessel_regions import masks,path_coordinates,section_definitions,JUMP
from flow_solver_support.wss_case import recover,material
from flow_solver_support.flow_parser import parse_solver_log

def csvout(p,rows):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if not rows:return
 with p.open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
def quant(v,a,qs=(.05,.5,.95)):
 order=np.argsort(v);return np.interp(qs,(np.cumsum(a[order])-.5*a[order])/a.sum(),v[order])
def stats(w,a,c):
 q=quant(w,a);mini=int(np.argmin(w));maxi=int(np.argmax(w));row=dict(facets=len(w),area_um2=float(a.sum()*1e12),mean_Pa=float(np.average(w,weights=a)),p05_Pa=float(q[0]),p50_Pa=float(q[1]),p95_Pa=float(q[2]),min_Pa=float(w[mini]),max_Pa=float(w[maxi]))
 for lab,xyz in [('min',c[mini]),('max',c[maxi])]:
  for j,axis in enumerate('xyz'):row[lab+'_'+axis+'_um']=float(xyz[j]*1e6)
 for lab,mask in [('below1',w<1),('below2',w<2),('above30',w>30)]:
  row[lab+'_area_pct']=float(100*a[mask].sum()/a.sum())
  centroid=np.average(c[mask],axis=0,weights=a[mask])*1e6 if mask.any() else [None]*3
  for j,axis in enumerate('xyz'):row[lab+'_centroid_'+axis+'_um']=float(centroid[j]) if mask.any() else ''
 return row

def analyze(case):
 case=Path(case).resolve();out=case/'reports';out.mkdir(exist_ok=True);mu=material(case)['mu_Pa_s'];a=np.load(case/'SV_MESH/mesh_arrays.npz');f=np.load(case/'frozen_flow/flow_arrays_si.npz');x=a['points_m'];t=a['tetra'];b=a['boundary_triangles'];tags=a['facet_tags'];u=f['velocity_m_s'];p=f['pressure_pa'];assert np.array_equal(x,f['points_m']) and np.array_equal(t,f['tetra'])
 for key in ['velocity_m_s','pressure_pa']:assert np.isfinite(f[key]).all()
 cmd=[sys.executable,'-B',str(C/'rotate_visualization/prepare_surface_data.py'),'--case',str(case),'--output',str(case/'wss')]
 with (out/'production_wss.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
 surf,d=recover(x,t,b,tags,u,p,mu);w=d['magnitude'];area=d['area'];c=d['centers'];normal=d['normal'];own=d['owners'];wall=b[tags==1];reg=masks(c)
 vol=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])/6;h=3*vol[own]/area;q=a['min_sicn'];policy=json.loads((case/'policy.json').read_text());meshname=policy.get('mesh_reference',case.name)
 rows=[]
 for name,sel in reg.items():
  if not sel.any():raise ValueError('Empty fixed physical region '+name)
  row=dict(case=case.name,mesh=meshname,region=name,statistical_weight='wall_triangle_area',mask_definition='fixed_physical_coordinates',reference=('H0_geometry_only_O2_changed_from_'+policy['comparison_baseline_case'] if policy.get('perturbation') else 'H0_geometry_H0_BC'),**stats(w[sel],area[sel],c[sel]),nearwall_mean_um=float(np.average(h[sel],weights=area[sel])*1e6),nearwall_p05_um=float(quant(h[sel],area[sel])[0]*1e6),nearwall_p50_um=float(quant(h[sel],area[sel])[1]*1e6),nearwall_p95_um=float(quant(h[sel],area[sel])[2]*1e6),owner_minSICN_min=float(q[own[sel]].min()),owner_minSICN_p05=float(np.quantile(q[own[sel]],.05)),owner_minSICN_median=float(np.median(q[own[sel]])))
  rows.append(row)
 csvout(out/'region_wss.csv',rows)
 baseline=V/'stage3/vessel_baseline/wss/data/wall_wss_si.vtp';fixed=[];overlap=[]
 if baseline.exists():
  base=pv.read(baseline);index,closest=base.find_closest_cell(c,return_closest_point=True);dist=np.linalg.norm(c-closest,axis=1);bw=np.asarray(base.cell_data['WSS_raw_Pa'])[index]
  for name,sel in reg.items():
   for maskname,mask in [('baseline_below1',bw<1),('baseline_below2',bw<2),('baseline_above30',bw>30)]:
    m=sel&mask
    if m.any():fixed.append(dict(case=case.name,mesh=meshname,region=name,statistical_weight='current_triangle_area',mask_definition=maskname+'_nearest_baseline_triangle',reference='vessel_baseline_frozen_mask',mapping_distance_max_um=float(dist[m].max()*1e6),mapping_distance_p95_um=float(np.quantile(dist[m],.95)*1e6),**stats(w[m],area[m],c[m])))
  csvout(out/'fixed_baseline_mask_wss.csv',fixed)
  for name,sel in reg.items():
   for threshold,old,new in [('below1',bw<1,w<1),('below2',bw<2,w<2),('above30',bw>30,w>30)]:
    inter=float(area[sel&old&new].sum());union=float(area[sel&(old|new)].sum());den=float(area[sel&old].sum()+area[sel&new].sum())
    overlap.append(dict(case=case.name,mesh=meshname,region=name,threshold=threshold,current_area_pct=float(100*area[sel&new].sum()/area[sel].sum()),mapped_baseline_area_pct=float(100*area[sel&old].sum()/area[sel].sum()),intersection_area_um2=inter*1e12,union_area_um2=union*1e12,Jaccard=inter/union if union else '',Dice=2*inter/den if den else '',statistical_weight='current_wall_triangle_area',reference='closest_baseline_facet_mask; empty union is undefined'))
  csvout(out/'threshold_overlap.csv',overlap)
 bx=x[b];av=.5*np.cross(bx[:,1]-bx[:,0],bx[:,2]-bx[:,0]);facearea=np.linalg.norm(av,axis=1);fl=np.einsum('ij,ij->i',av,u[b].mean(axis=1));flows={ROLES[int(tag)]:float(fl[tags==tag].sum()) for tag in np.unique(tags)};qin=-flows['INLET'];boundaryrows=[]
 for tag,role in ROLES.items():
  sel=tags==tag
  boundaryrows.append(dict(case=case.name,mesh=meshname,boundary=role,Q_outward_m3_s=flows[role],Q_uL_min=flows[role]*6e10,flow_relative_to_inlet_pct=100*flows[role]/qin,pressure_area_mean_Pa=float(np.average(p[b[sel]].mean(axis=1),weights=facearea[sel])),statistical_weight='exact_P1_triangle_flux;area_pressure',reference='H0_boundary_roles'))
 csvout(out/'boundary_flows.csv',boundaryrows)
 # Pointwise recovered traction diagnostic is distinct from the weak FEM residual.
 caprows=[];boundary_owner=wss.boundary_owners(t,b);bn=av/facearea[:,None];gg=d['gradient'][boundary_owner];visc=mu*np.einsum('nij,nj->ni',gg+gg.swapaxes(1,2),bn);normalvisc=np.einsum('ij,ij->i',visc,bn);pressureface=p[b].mean(axis=1);tree=ET.parse(case/'run/solver.xml')
 for tag in [2,3,5]:
  role=ROLES[tag];sel=tags==tag;prescribed=float(tree.find('.//Add_BC[@name="%s"]/Value'%role).text);res=normalvisc[sel]-pressureface[sel]+prescribed
  caprows.append(dict(case=case.name,mesh=meshname,boundary=role,prescribed_pressure_Pa=prescribed,pressure_area_mean_Pa=float(np.average(pressureface[sel],weights=facearea[sel])),normal_viscous_stress_area_mean_Pa=float(np.average(normalvisc[sel],weights=facearea[sel])),normal_traction_residual_mean_Pa=float(np.average(res,weights=facearea[sel])),normal_traction_residual_RMS_Pa=float(np.sqrt(np.average(res**2,weights=facearea[sel]))),statistical_weight='cap_triangle_area',reference='P1 recovered sigma*n minus prescribed -P*n; not assembled weak residual'))
 csvout(out/'cap_traction.csv',caprows)
 # Exact P1 interpolation onto planar intersection polygons, triangulated for integration.
 grid=pv.UnstructuredGrid(np.column_stack([np.full(len(t),4),t]).ravel(),np.full(len(t),pv.CellType.TETRA,np.uint8),x);grid['Velocity']=u;grid['Pressure']=p
 sections=[]
 for spec in section_definitions():
  n=np.array(spec['normal']);origin=np.array(spec['origin_m']);sl=grid.slice(normal=n,origin=origin,generate_triangles=True);tri=sl.faces.reshape(-1,4)[:,1:];xyz=sl.points[tri];cent=xyz.mean(axis=1);sel=np.linalg.norm(cent-origin,axis=1)<spec['radius_limit_m'];tri=tri[sel];xyz=xyz[sel];aa=.5*np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1);qsec=float(np.dot(aa,np.asarray(sl['Velocity'])[tri].mean(axis=1)@n));ref=(flows['OUTLET_01']+flows['OUTLET_03']) if spec['compare_boundary']=='O1_plus_O3' else flows['INLET' if spec['compare_boundary']=='INLET' else 'OUTLET_0'+spec['compare_boundary'][1]];ref*=spec['expected_sign']
  connected=sl.extract_cells(np.flatnonzero(sel)).connectivity()
  components=len(np.unique(connected.cell_data['RegionId']))
  radiusmax=float(np.linalg.norm(xyz-origin,axis=-1).max())
  assert components==1 and radiusmax<spec['radius_limit_m'],'Section must be one lumen entirely inside the selection disk'
  sections.append(dict(case=case.name,mesh=meshname,section=spec['name'],**{k:json.dumps(v) if isinstance(v,list) else v for k,v in spec.items() if k!='name'},area_um2=float(aa.sum()*1e12),selected_connected_components=components,selected_max_vertex_radius_um=radiusmax*1e6,selection_disk_does_not_cut_lumen=True,Q_m3_s=qsec,reference_boundary_Q_m3_s=ref,Q_difference_pct_of_inlet=100*(qsec-ref)/qin,pressure_area_mean_Pa=float(np.average(np.asarray(sl['Pressure'])[tri].mean(axis=1),weights=aa)),statistical_weight='exact_P1_planar_triangle_flux',reference='fixed_physical_plane_and_radius'))
 csvout(out/'internal_sections.csv',sections)
 edges=np.sort(wall[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1);_,inv,count=np.unique(edges,axis=0,return_inverse=True,return_counts=True);order=np.argsort(inv,kind='stable');starts=np.r_[0,np.cumsum(count)[:-1]];ii=starts[count==2];left=order[ii]//3;right=order[ii+1]//3;mid=(c[left]+c[right])/2;jump=np.abs(w[left]-w[right]);angle=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',normal[left],normal[right]),-1,1)));dist=np.linalg.norm(mid-JUMP,axis=1);adj=[]
 for name,sel in [('all_wall',np.ones(len(left),bool)),('J1_r5um',np.linalg.norm(mid-np.array([92,49,111])*1e-6,axis=1)<5e-6),('known_jump_r1um',dist<1e-6)]:
  vals=jump[sel];kk=np.flatnonzero(sel)[np.argmax(vals)];adj.append(dict(case=case.name,mesh=meshname,region=name,pairs=int(sel.sum()),jump_mean_Pa=float(vals.mean()),jump_p50_Pa=float(np.quantile(vals,.5)),jump_p95_Pa=float(np.quantile(vals,.95)),jump_max_Pa=float(vals.max()),max_mid_x_um=float(mid[kk,0]*1e6),max_mid_y_um=float(mid[kk,1]*1e6),max_mid_z_um=float(mid[kk,2]*1e6),statistical_weight='adjacent_pair_count',reference='shared_wall_edges'))
 csvout(out/'adjacent_jump_summary.csv',adj);detail=[]
 for k in sorted(set(np.argsort(dist)[:12].tolist()+np.argsort(jump*(dist<1e-6))[-12:].tolist())):
  l=int(left[k]);r=int(right[k]);detail.append(dict(case=case.name,mesh=meshname,facet_a=int(d['wall_ids'][l]),facet_b=int(d['wall_ids'][r]),parent_a=int(own[l]),parent_b=int(own[r]),x_um=float(mid[k,0]*1e6),y_um=float(mid[k,1]*1e6),z_um=float(mid[k,2]*1e6),distance_to_original_jump_um=float(dist[k]*1e6),wss_a_Pa=float(w[l]),wss_b_Pa=float(w[r]),jump_Pa=float(jump[k]),normal_angle_deg=float(angle[k]),sicn_a=float(q[own[l]]),sicn_b=float(q[own[r]]),reference='physical_location_neighborhood_not_old_ID'))
 csvout(out/'known_jump_pairs.csv',detail)
 sbest,pmask,total=path_coordinates(c);profile=[]
 for lo in np.arange(0,total*1e6,1):
  sel=pmask&(sbest*1e6>=lo)&(sbest*1e6<lo+1)
  if sel.any():profile.append(dict(case=case.name,mesh=meshname,region='O2_fixed_path',s_mid_um=float(lo+.5),statistical_weight='wall_triangle_area',reference='audit_fixed_path_1um_bins',**stats(w[sel],area[sel],c[sel])))
 csvout(out/'outlet2_profile.csv',profile)
 gradient=d['gradient'];div=np.trace(gradient,axis1=1,axis2=2);quality=dict(epsilon_mass=abs(sum(flows.values()))/qin,epsilon_Q=abs(qin-policy['Q_target_m3_s'])/policy['Q_target_m3_s'],wall_max_speed_m_s=float(np.linalg.norm(u[np.unique(wall)],axis=1).max()),volume_RMS_divergence_s_inv=float(np.sqrt(np.average(div**2,weights=vol))),RMS_divergence_over_RMS_gradient=float(np.sqrt(np.dot(vol,div**2)/np.dot(vol,np.sum(gradient**2,axis=(1,2))))))
 wallflag=np.zeros(len(x),bool);wallflag[np.unique(wall)]=True;allwall=np.all(wallflag[t],axis=1);quality.update(tetra_all_four_nodes_on_wall=int(allwall.sum()),wall_facets_with_all_wall_owner=int(allwall[own].sum()),exact_zero_WSS_facets=int(np.count_nonzero(w==0)),exact_zero_WSS_area_pct=float(100*area[w==0].sum()/area.sum()))
 if allwall.any():
  bad=[]
  for k in np.flatnonzero(allwall):
   xyz=x[t[k]].mean(axis=0)*1e6;bad.append(dict(case=case.name,cell_zero_based=int(k),x_um=float(xyz[0]),y_um=float(xyz[1]),z_um=float(xyz[2]),minSICN=float(q[k]),volume_um3=float(vol[k]*1e18),reason='all four P1 velocity nodes lie on no-slip wall'))
  csvout(out/'all_wall_tetra.csv',bad)
 logpath=case/'run/solver.log'
 if logpath.exists():
  log=parse_solver_log(logpath.read_text(errors='replace'),policy['dt_s']);solves=log['linear_solves'];rat=[]
  for s in solves:
   ms=s.get('petsc_monitor',[])
   if ms:rat.append(ms[-1]['true_residual_norm']/max(1e-24,1e-10*ms[0]['true_residual_norm']))
  quality.update(parsed_solves=len(solves),unparsed_rows=len(log['unparsed_rows']),failed_linear_solves=log['failed_linear_solves'],nonlinear_failure_messages=log['nonlinear_failure_messages'],recovered_attempts=log['recovered_attempts'],last_nonlinear_relative=solves[-1]['nonlinear_Ri_over_R0'],max_true_residual_over_criterion=max(rat) if rat else None)
 quality['accepted_final_and_log_checks']=bool(quality.get('parsed_solves',0) and quality['unparsed_rows']==0 and not quality['failed_linear_solves'] and not quality['nonlinear_failure_messages'] and quality['last_nonlinear_relative']<1e-10 and quality['max_true_residual_over_criterion']<1.01 and quality['epsilon_mass']<=1e-6 and quality['epsilon_Q']<=1e-6 and quality['wall_max_speed_m_s']<=1e-12)
 dump(out/'flow_quality.json',quality);print(json.dumps(dict(case=case.name,regions=len(rows),quality=quality),indent=2),flush=True)
 assert quality['accepted_final_and_log_checks'],'Actual final field/log gate failed; retain diagnostics, do not accept comparison as converged CFD'
 return rows
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);args=p.parse_args();analyze(args.case)
