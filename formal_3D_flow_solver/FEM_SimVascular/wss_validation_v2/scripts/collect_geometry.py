"""Source-surface differences kept separate from numerical mesh refinement."""
from case_common import *
from analyze_vessel import csvout
rows=[]
for case in [V/'stage3'/name for name in ['vessel_baseline','vessel_medium','vessel_fine']]:
 q=case/'reports/geometry_qc.json';old=case/'reports/formal_geometry_qc.json'
 if q.exists():
  a=json.loads(q.read_text());v0=a['source_volume_m3'];v1=a['actual_volume_m3'];d=a['wall_bidirectional_distance_m'];ports=a['ports'];h=json.loads((case/'mesh_request.json').read_text())['global_edge_size_m']
 elif old.exists():
  a=json.loads(old.read_text());v0=a['source_enclosed_volume_m3'];v1=a['actual_enclosed_volume_m3'];d=dict(max=a['wall']['max_m'],p95=a['wall']['P95_m'],rms=a['wall']['RMS_m']);ports={k:dict(source=v['original'],actual=v['actual'],area_relative_change=v['actual']['area_m2']/v['original']['area_m2']-1,centroid_shift_m=v['centroid_shift_m'],normal_dot=v['normal_dot']) for k,v in a['ports'].items()};h=json.loads((C/'configs/mesh_policy.json').read_text())['h_ref_m']
 else:continue
 for role,pp in ports.items():
  rows.append(dict(case=case.name,mesh=case.name,CFD_status=('SKIPPED_BY_USER' if (case/'reports/user_cancellation.json').exists() else 'QUALIFIED_REUSED_CFD' if (case/'reports/baseline_reuse.json').exists() else json.loads((case/'reports/execution.json').read_text()).get('status','NOT_RUN') if (case/'reports/execution.json').exists() else 'NOT_YET_COMPLETED'),region=role,global_h_um=h*1e6,source_volume_um3=v0*1e18,actual_volume_um3=v1*1e18,volume_signed_difference_pct=100*(v1/v0-1),source_port_area_um2=pp['source']['area_m2']*1e12,actual_port_area_um2=pp['actual']['area_m2']*1e12,port_area_signed_difference_pct=100*pp['area_relative_change'],source_equivalent_radius_um=pp['source']['equivalent_radius_m']*1e6,actual_equivalent_radius_um=pp['actual']['equivalent_radius_m']*1e6,port_centroid_shift_um=pp['centroid_shift_m']*1e6,port_normal_dot=pp['normal_dot'],wall_bidirectional_sample_max_um=d['max']*1e6,wall_bidirectional_sample_p95_um=d['p95']*1e6,wall_bidirectional_sample_RMS_um=d['rms']*1e6,statistical_weight='triangle_area_for_ports;all_vertices_and_triangle_centers_for_distance',reference='same_source_exterior_surface_before_TetGen'))
if rows:csvout(V/'data/vessel_geometry_comparison.csv',rows)
print(len(rows),'port geometry comparison rows')
