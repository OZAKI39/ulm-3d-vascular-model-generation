"""Full signed-traction and magnitude responses on the actually solved fixed mesh.

Run only after stage four CFD and production WSS have passed independent checks.
"""
import numpy as np
import pyvista as pv
from case_common import *
from vessel_regions import masks
from analyze_vessel import csvout

group=json.loads((V/'stage4/boundary_group.json').read_text())
base=Path(group['baseline_path'])
cases=[base,V/'stage4/O2_minus1pct',V/'stage4/O2_plus1pct']
surfaces=[]
for case in cases:
    assert json.loads((case/'reports/execution.json').read_text())['status']=='PASS'
    assert json.loads((case/'reports/flow_quality.json').read_text())['accepted_final_and_log_checks']
    surfaces.append(pv.read(case/'wss/data/wall_wss_si.vtp'))
s0,sm,sp=surfaces
for s in [sm,sp]:
    assert np.array_equal(s.points,s0.points) and np.array_equal(s.faces,s0.faces)
    for key in ['Global_boundary_facet_zero_based','Area_m2','Outward_normal','Parent_tetra_zero_based']:
        assert np.array_equal(s.cell_data[key],s0.cell_data[key])
area=np.asarray(s0.cell_data['Area_m2']);cent=s0.cell_centers().points
t0,tm,tp=[np.asarray(s.cell_data['Tangential_viscous_traction_Pa']) for s in surfaces]
w0,wm,wp=[np.asarray(s.cell_data['WSS_raw_Pa']) for s in surfaces]
rows=[]
for region,sel in masks(cent).items():
    aa=area[sel]
    def rms(z):
        v=np.sum(z[sel]**2,axis=1) if z.ndim==2 else z[sel]**2
        return float(np.sqrt(np.average(v,weights=aa)))
    reference=rms(w0);vector_signal=max(rms(tm-t0),rms(tp-t0));magnitude_signal=max(rms(wm-w0),rms(wp-w0))
    for label,t,w in [('minus1pct',tm,wm),('plus1pct',tp,wp)]:
        rows.append(dict(baseline_case=base.name,case='O2_'+label,mesh=group['mesh'],region=region,delta_O2_pressure_Pa=(-1 if label=='minus1pct' else 1)*group['delta_Pa'],reference_WSS_area_RMS_Pa=reference,traction_vector_response_area_RMS_Pa=rms(t-t0),traction_vector_response_relative_pct=100*rms(t-t0)/reference,WSS_magnitude_response_area_RMS_Pa=rms(w-w0),WSS_magnitude_response_relative_pct=100*rms(w-w0)/reference,vector_even_response_area_RMS_Pa=rms(tp+tm-2*t0),vector_symmetry_defect_over_larger_response=rms(tp+tm-2*t0)/vector_signal if vector_signal else '',magnitude_even_response_area_RMS_Pa=rms(wp+wm-2*w0),magnitude_symmetry_defect_over_larger_response=rms(wp+wm-2*w0)/magnitude_signal if magnitude_signal else '',statistical_weight='same_medium_wall_triangle_area',reference='same_mesh_actual_CFD_baseline; signed vectors separated from magnitudes'))
csvout(V/'data/boundary_field_responses.csv',rows)
print('Actual fixed-grid full-field response rows:',len(rows))
