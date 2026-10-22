"""Identical-mesh backend/time-step metric changes, separate from spatial sequence."""
from pathlib import Path
import csv,json
V=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((V/'data/vessel_region_wss.csv').open()));look={(r['case'],r['region']):r for r in rows};out=[]
for case,kind in [('vessel_baseline_cpu_mpi8','execution_backend'),('vessel_baseline_gpu_mpi1_halfdt','time_step')]:
 for (name,region),b in look.items():
  if name!='vessel_baseline' or (case,region) not in look:continue
  c=look[(case,region)]
  for metric in ['mean_Pa','p05_Pa','p50_Pa','p95_Pa','min_Pa','max_Pa','below1_area_pct','below2_area_pct','above30_area_pct']:
   a,z=float(b[metric]),float(c[metric]);out.append(dict(reference_case='vessel_baseline',comparison_case=case,comparison_type=kind,mesh='vessel_baseline',region=region,metric=metric,unit='percentage_points' if 'area_pct' in metric else 'Pa',reference_value=a,comparison_value=z,signed_absolute_change=z-a,signed_relative_change_pct=100*(z/a-1) if a else '',statistical_weight='same_wall_triangle_area',reference='identical mesh; original H0'))
with (V/'data/vessel_control_metric_changes.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
print('Saved',len(out),'control comparisons, independent of the main three-grid sequence.')
