"""Compact, measured two-mesh comparisons; no extrapolation or new CFD."""
from pathlib import Path
import csv,json,math
V=Path(__file__).resolve().parents[1]
def read(p):return list(csv.DictReader(Path(p).open()))
def write(name,rows):
 with (V/'data'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
names=['vessel_baseline','vessel_medium']
cases=[V/'stage3'/n for n in names]
for c in cases:assert json.loads((c/'reports/flow_quality.json').read_text())['accepted_final_and_log_checks']
rr=[{r['region']:r for r in read(c/'reports/region_wss.csv')} for c in cases]
spatial=[]
for region in rr[0]:
 a,b=[r[region] for r in rr]
 for label in ['min','max','below1_centroid','below2_centroid','above30_centroid']:
  keys=[label+'_'+ax+'_um' for ax in 'xyz'];x=[a[k] for k in keys];y=[b[k] for k in keys]
  valid=all(z!='' for z in x+y)
  distance=math.sqrt(sum((float(z)-float(w))**2 for z,w in zip(x,y))) if valid else ''
  spatial.append(dict(region=region,position_type=label,coordinate_unit='um',baseline_case=names[0],comparison_case=names[1],
    **{f'baseline_{ax}_um':z for ax,z in zip('xyz',x)},**{f'medium_{ax}_um':z for ax,z in zip('xyz',y)},
    displacement_um=distance,reference='fixed_physical_region; centroids weighted by each mesh own threshold area',
    limitation='Empty threshold has no centroid. Centroid displacement alone does not measure patch shape or validate a material-point correspondence.'))
write('vessel_spatial_changes.csv',spatial)
ss=[{r['section']:r for r in read(c/'reports/internal_sections.csv')} for c in cases];rows=[]
for section in ss[0]:
 a,b=[float(z[section]['Q_difference_pct_of_inlet']) for z in ss]
 rows.append(dict(section=section,baseline_case=names[0],comparison_case=names[1],baseline_signed_deviation_pct_of_inlet=a,medium_signed_deviation_pct_of_inlet=b,change_percentage_points=b-a,absolute_deviation_reduction_pct=100*(1-abs(b)/abs(a)) if a else '',statistical_weight='exact_P1_planar_flux',reference='same physical plane; corresponding cap flow or O1+O3 flow'))
write('vessel_internal_flow_comparison.csv',rows)
ff=[{r['boundary']:r for r in read(c/'reports/boundary_flows.csv')} for c in cases];flows=[]
for boundary in ['OUTLET_01','OUTLET_02','OUTLET_03']:
 a,b=[r[boundary] for r in ff];q0,q1=[float(r['Q_outward_m3_s']) for r in [a,b]];f0,f1=[float(r['flow_relative_to_inlet_pct']) for r in [a,b]]
 flows.append(dict(boundary=boundary,baseline_case=names[0],comparison_case=names[1],baseline_Q_m3_s=q0,medium_Q_m3_s=q1,baseline_fraction_pct=f0,medium_fraction_pct=f1,change_percentage_points=f1-f0,flow_relative_change_pct=100*(q1/q0-1),statistical_weight='exact_P1_boundary_triangle_flux',reference='unchanged inlet and outlet BCs; two meshes'))
write('vessel_flow_split_comparison.csv',flows)
aa=[{r['region']:r for r in read(c/'reports/adjacent_jump_summary.csv')} for c in cases];jumps=[]
for region in aa[0]:
 for key in ['jump_mean_Pa','jump_p50_Pa','jump_p95_Pa','jump_max_Pa']:
  a,b=[float(r[region][key]) for r in aa]
  jumps.append(dict(region=region,metric=key,unit='Pa',baseline_case=names[0],comparison_case=names[1],baseline_value=a,medium_value=b,signed_absolute_change_Pa=b-a,relative_change_pct=100*(b/a-1),statistical_weight='shared_edge_pair_count',reference='same physical region, not same facet IDs'))
write('vessel_local_jump_comparison.csv',jumps)
write('validation_stage_status.csv',[
 dict(stage='1',object='production WSS regression and discriminating tests',status='COMPLETED_ACTUAL_CHECKS',reason='Production entry executed; raw/display maximum difference 0 Pa'),
 dict(stage='2',object='three-mesh pipe CFD and fine pipe half-dt',status='COMPLETED_ACTUAL_CFD',reason='Kept completed results; no repetition after scope restriction'),
 dict(stage='3',object='original and medium vascular CFD',status='COMPLETED_TWO_MESH_SENSITIVITY',reason='Identity-checked original CFD and actual accepted medium solve; no mesh-independence claim'),
 dict(stage='3',object='fine vascular CFD',status='SKIPPED_BY_USER',reason='Cancelled for time cost before initialization/CFD; mesh and configuration retained'),
 dict(stage='4',object='outlet pressure sensitivity CFD',status='SKIPPED_BY_USER',reason='因用户限制时间成本，未执行; cases not generated, queued or launched')])
tables=[]
tables+=['| 区域 | 均值Pa：原→中 | P5 Pa：原→中 | P50 Pa：原→中 | P95 Pa：原→中 | 均值变化% |','|---|---:|---:|---:|---:|---:|']
for region in rr[0]:
 a,b=[r[region] for r in rr];values=['%.4f → %.4f'%(float(a[k]),float(b[k])) for k in ['mean_Pa','p05_Pa','p50_Pa','p95_Pa']]
 tables.append('| '+region+' | '+' | '.join(values)+' | %.3f |'%(100*(float(b['mean_Pa'])/float(a['mean_Pa'])-1)))
tables+=['','| 区域 | <1Pa面积%：原→中 | <2Pa面积%：原→中 | >30Pa面积%：原→中 |','|---|---:|---:|---:|']
for region in rr[0]:
 a,b=[r[region] for r in rr];tables.append('| '+region+' | '+' | '.join('%.4f → %.4f'%(float(a[k]),float(b[k])) for k in ['below1_area_pct','below2_area_pct','above30_area_pct'])+' |')
tables+=['','| 固定截面 | 原网格偏差%Qin | 中档偏差%Qin | 绝对偏差降低% |','|---|---:|---:|---:|']
for r in rows:tables.append('| %s | %.5f | %.5f | %.2f |'%(r['section'],r['baseline_signed_deviation_pct_of_inlet'],r['medium_signed_deviation_pct_of_inlet'],r['absolute_deviation_reduction_pct']))
tables+=['','| 出口 | 原流量m³/s | 中流量m³/s | 原分流% | 中分流% | 差：百分点 |','|---|---:|---:|---:|---:|---:|']
for r in flows:tables.append('| %s | %.8e | %.8e | %.6f | %.6f | %+.6f |'%(r['boundary'],r['baseline_Q_m3_s'],r['medium_Q_m3_s'],r['baseline_fraction_pct'],r['medium_fraction_pct'],r['change_percentage_points']))
(V/'data/two_mesh_report_tables.md').write_text('\n'.join(tables)+'\n')
print('Saved two-mesh spatial, local-jump, internal-flux, flow-split, stage-status CSVs and report tables.')
