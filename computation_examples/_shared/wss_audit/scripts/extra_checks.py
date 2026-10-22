"""Additional read-only BC, local mesh, pressure sign, and geometry evidence."""
import sys, json, csv
from pathlib import Path
import numpy as np
sys.dont_write_bytecode=True
from audit import A,B,C,H,O,V,N,ND,prod,dump,csvout,sha
sys.path.insert(0,str(N))
from network_1d0d.network_solver import solve_network,operating_point_scale
from network_1d0d.hydraulic_resistance import linear_radius_resistance,ROI_TARGET_Q_M3_S
from network_1d0d.boundary_conditions import shift_pressure_gauge,cap_pressures

def main():
    g=np.load(A/'inputs/network/analysis_A_H0_graph_si.npz');x=g['xyz_m'];r=g['radius_m'];e=g['edges'];src=int(g['source_index'])
    resistance=linear_radius_resistance(np.linalg.norm(x[e[:,1]]-x[e[:,0]],axis=1),r[e[:,0]],r[e[:,1]])
    bc={int(i):0. for i in g['terminal_indices']};bc[src]=1.
    sol=solve_network(len(x),e,resistance,bc);inlet=int(g['port_indices'][0])
    ei=np.flatnonzero(g['roi_internal_edge_mask']&np.any(e==inlet,axis=1));assert len(ei)==1
    qin=sol.flow[ei[0]]*(1 if e[ei[0],0]==inlet else -1);scale=operating_point_scale(qin,ROI_TARGET_Q_M3_S)
    saved=np.load(ND/'H0_solution_si.npz');p_error=abs(sol.pressure*scale-saved['pressure_Pa']).max();q_error=abs(sol.flow*scale-saved['edge_Q_m3s']).max()
    transfer=json.loads((A/'inputs/network/roi_fixed_pressure_bc_H0.json').read_text())
    rows=transfer['ports'];raw=np.array([cap_pressures(a['pressure_realcut_Pa'],a['signed_Q_m3s'],a['R_extension_Pa_s_m3']) for a in rows]);shift=shift_pressure_gauge(raw)
    record=dict(source_node_id=int(g['ids'][src]),terminals=len(g['terminal_indices']),nodes=len(x),edges=len(e),unit_source_pressure_Pa=1.,unit_roi_Q_m3_s=qin,pressure_scale=scale,actual_source_pressure_Pa=scale,
        pressure_reproduction_max_error_Pa=p_error,flow_reproduction_max_error_m3_s=q_error,network_solve_audit=sol.audit,
        cap_raw_Pa=raw,cap_shifted_Pa=shift,pressure_transfer_max_error_Pa=abs(shift-np.array([a['pressure_cap_shifted_Pa'] for a in rows])).max(),
        pressure_pairwise_shift_error_Pa=abs((raw[:,None]-raw)-(shift[:,None]-shift)).max(),note='Network solve and pressure transfer actually rerun; no 3D solver run.')
    dump('data/network_reproduction.json',record)
    d=np.load(A/'inputs/H0/flow_arrays_si.npz');m=np.load(A/'inputs/mesh_arrays.npz');pa=np.load(A/'data/plot_arrays.npz');xyz=d['points_m'];t=d['tetra'];u=d['velocity_m_s'];p=d['pressure_pa'];b=d['boundary_triangles'];tags=d['facet_tags'];mu=.00345312
    grad=prod.p1_gradients(xyz,t,u);allown=prod.boundary_owners(t,b);cent,area,normal=prod.wall_geometry(xyz,t,b,allown)
    summary=[]
    for tag,role,pbc in [(4,'INLET',None),(3,'O1',585.2864332598318),(5,'O2',2932.0152710782013),(2,'O3',0.)]:
        s=tags==tag;uv=u[b[s]].mean(axis=1);un=(uv*normal[s]).sum(1);normalvis=mu*np.einsum('ni,nij,nj->n',normal[s],grad[allown[s]]+grad[allown[s]].swapaxes(1,2),normal[s]);avp=np.average(p[b[s]].mean(1),weights=area[s])
        summary.append(dict(region=role,cap_area_um2=area[s].sum()*1e12,cap_radius_um=np.sqrt(area[s].sum()/np.pi)*1e6,cap_mean_pressure_Pa=avp,cap_mean_normal_viscous_Pa=np.average(normalvis,weights=area[s]),prescribed_pressure_Pa=pbc,
            reverse_flux_facets=int((un>0).sum()) if role=='INLET' else int((un<0).sum()),max_normal_speed_m_s=abs(un).max()))
    csvout('data/cap_traction_check.csv',summary)
    w=pa['wss'];centers=pa['centers_um'];owns=pa['own'];wallids=np.flatnonzero(tags==1);iq=int(np.argmin(m['min_sicn']));ext=[]
    for label,ii in [('global_min_WSS',np.argmin(w)),('global_max_WSS',np.argmax(w))]:
        ext.append(dict(label=label,facet_zero_based=int(wallids[ii]),tetra_zero_based=int(owns[ii]),xyz_um=centers[ii],wss_Pa=w[ii],minSICN=m['min_sicn'][owns[ii]],altitude_um=pa['h_um'][ii]))
    qi=xyz[t[iq]].mean(0)*1e6;mind=float(np.min(np.linalg.norm(centers[w<2]-qi,axis=1)))
    ext.append(dict(label='minimum_mesh_SICN',tetra_zero_based=iq,xyz_um=qi,minSICN=m['min_sicn'][iq],distance_to_first_junction_um=np.linalg.norm(qi-pa['first_junction_um']),nearest_low_WSS_center_distance_um=mind,is_wall_owner=bool(np.any(owns==iq))))
    # Exact 1-hop vicinity of minimum WSS, plus curvature/normal consistency is handled in main.
    dump('data/critical_locations.json',ext)
    # WSS / Poiseuille comparison is a scale check in artificial straight extensions, not a proof.
    regions=list(csv.DictReader((A/'data/region_summary.csv').open()));bound=list(csv.DictReader((A/'data/boundary_values.csv').open()));comp=[]
    for i in [1,2,3]:
        q=float(next(a['Q_outward_m3_s'] for a in bound if a['case']=='H0' and a['boundary']==f'OUTLET_0{i}'))
        radius=next(a['cap_radius_um'] for a in summary if a['region']==f'O{i}')*1e-6
        expected=4*mu*q/(np.pi*radius**3);actual=float(next(a['mean_Pa'] for a in regions if a['case']=='H0' and a['region']==f'O{i}_extension'))
        comp.append(dict(region=f'O{i}_extension',case='H0',Q_m3_s=q,cap_equivalent_R_um=radius*1e6,Poiseuille_scale_Pa=expected,computed_area_mean_Pa=actual,relative_difference_pct=100*(actual/expected-1),limitation='Cap equivalent radius; noncircular polygon and entrance/end effects; not exact analytic solution'))
    csvout('data/extension_scale_check.csv',comp)
    # Record local vs remote solver/config/field identity; remote command output preserved separately.
    print(json.dumps(record,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item(),indent=2)); print(ext)

if __name__=='__main__':main()
