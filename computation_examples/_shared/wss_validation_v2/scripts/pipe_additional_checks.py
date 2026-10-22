"""Explain observed P1 recovery bias and report internal mass/gradient evidence."""
import csv,json
import numpy as np
from pathlib import Path
from case_common import *
from analyze_vessel import csvout
geom=[];sections=[];divs=[]
for name in ['pipe_nr4','pipe_nr8','pipe_nr16_cpu_mpi8','pipe_nr16_cpu_mpi8_halfdt']:
 c=V/'stage2'/name
 if not (c/'reports/pipe_validation.csv').exists():continue
 p=json.loads((c/'policy.json').read_text());r=next(r for r in csv.DictReader((c/'reports/pipe_validation.csv').open()) if r['field_type']=='analytic_nodal');R=p['R_m'];dr=R/p['radial_intervals'];theta=np.pi/p['circumferential_intervals'];tau=4*p['mu_Pa_s']*p['Umean_m_s']/R;expected=tau*(1-dr/(2*R))/np.cos(theta);actual=float(r['wss_area_mean_Pa']);m=np.load(c/'SV_MESH/mesh_arrays.npz');x=m['points_m'];t=m['tetra'];b=m['boundary_triangles'];tags=m['facet_tags'];own=m['adjacent_tetra'][tags==1];wall=b[tags==1];area=np.linalg.norm(np.cross(x[wall[:,1]]-x[wall[:,0]],x[wall[:,2]]-x[wall[:,0]]),axis=1)/2;vol=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])/6;h=3*vol[own]/area
 geom.append(dict(case=name,region='pipe_wall',method='closed_form_P1_slope_on_this_ring_extrusion_vs_production_analytic_node_result',statistical_weight='wall_triangle_area',reference_continuum_tau_Pa=tau,expected_P1_tau_Pa=float(expected),production_analytic_nodal_tau_Pa=actual,mean_abs_difference_Pa=abs(actual-expected),normal_height_expected_um=dr*np.cos(theta)*1e6,normal_height_actual_mean_um=float(np.average(h,weights=area)*1e6),normal_height_max_error_um=float(np.max(abs(h-dr*np.cos(theta)))*1e6),normal_height_over_R=float(np.average(h,weights=area)/R)))
 ss=list(csv.DictReader((c/'reports/sections.csv').open()));qin=float(ss[0]['Q_m3_s'])
 for row in ss:row.update(dt_s=p['dt_s'],Q_difference_pct_of_inlet=100*(float(row['Q_m3_s'])/qin-1),Q_difference_pct_of_continuum=100*(float(row['Q_m3_s'])/float(row['Q_theory_m3_s'])-1),reference='same_case_inlet_and_continuous_Poiseuille_Q');sections.append(row)
 f=np.load(c/'frozen_flow/flow_arrays_si.npz');G=wss.p1_gradients(x,t,f['velocity_m_s']);dv=np.trace(G,axis1=1,axis2=2);cent=x[t].mean(axis=1);L=p['L_m']
 for region,sel in [('all_volume',np.ones(len(t),bool)),('central_z_0.3L_to_0.7L',(cent[:,2]>=.3*L)&(cent[:,2]<=.7*L))]:divs.append(dict(case=name,region=region,dt_s=p['dt_s'],volume_um3=float(vol[sel].sum()*1e18),divergence_RMS_s_inv=float(np.sqrt(np.average(dv[sel]**2,weights=vol[sel]))),divergence_RMS_over_gradient_RMS=float(np.sqrt(np.dot(vol[sel],dv[sel]**2)/np.dot(vol[sel],np.sum(G[sel]**2,axis=(1,2))))),mean_signed_divergence_s_inv=float(np.average(dv[sel],weights=vol[sel])),statistical_weight='tetra_volume',reference='P1 velocity gradient; exact incompressible continuum div=0'))
for name,rr in [('pipe_P1_geometry_explanation.csv',geom),('pipe_internal_sections.csv',sections),('pipe_divergence.csv',divs)]:csvout(V/'data'/name,rr)
print('Additional pipe geometry/section/divergence evidence generated.')
