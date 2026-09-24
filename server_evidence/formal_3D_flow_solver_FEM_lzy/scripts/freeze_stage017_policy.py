#!/usr/bin/env python3
import json,sys,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,write_json,timestamp
r=ROOT/'reports/stage01_7';i=ROOT/'inputs/stage01_7';assert not (r/'acceptance_policy.json').exists()
b=json.loads((r/'baseline_recomputed.json').read_text());assert b['status']=='PASS'
g=json.loads((ROOT/'reports/stage01_6/acceptance_policy.json').read_text())['geometry']
policy={'stage':'1.7','timestamp':timestamp(),'search':{'search_factor':1.25,'grading_slope':1.0,'maximum_bracket_refinements':4,'maximum_H_over_R_eq':2.,'volume_feedback_refinement_factor':1.20,
 'initial_H':'Stage 1 medium bulk target read from metadata','formula':'h(d)=min(H,h_rim+grading_slope*d)','surface_selection':['minimum actual passing triangle count','larger H on ties'],
 'geometry_failure':'stop port search before quality decision','port_tie_break':'lexical port name'},
 'limits':{'maximum_surface_trials_per_port':8,'maximum_volume_iterations':5,'maximum_total_volume_meshes':5,
 'surface_trial_scope':'initial search plus all volume-feedback remeshes, per optimizer execution','volume_mesh_scope':'production plus independent determinism replay combined',
 'budget_interpretation':'Conservative interpretation pending optional clarification; no extra production or replay budget',
 'volume_wall_time_s':300,'volume_virtual_memory_gib':8},
 'geometry':g,'surface_quality':{'low_threshold':.1,'low_count_max':0,'P5_min':.45,'median_min':.70,'positive_finite_triangles':True},
 'surface_meshing':{'algorithm':6,'normalized_geometry_tolerance':1e-12,'random_seed':1,'normalization':'planar coordinates divided by current h_rim; original 3D rim reused exactly'},
 'volume_quality':{'low_threshold':.1,'cap_low_fraction':.30,'total_low_fraction':.50,'quantile_roundoff_tolerance':1e-12,'quantiles':['P1','P5','median'],'minimum_is_gate':False,
 'classification':'nearest exterior triangle center, same as Stage 1'},
 'cost':{'C_P2_max':1.35,'C_tetra_max':1.35,'cap_triangle_count_is_hard_gate':False},
 'volume_meshing':{'effective_gmsh_options':b['effective_gmsh_options'],'h_volume_m':b['h_volume_m'],'local_volume_refinement':False},
 'selection':['first feasible stops further production iterations','minimum P2 velocity proxy among actual feasible history','minimum tetra','minimum cap-adjacent low'],
 'determinism':{'surface':'repeat optimizer using same geometry/config/Gmsh; compare H, selected trial and canonical geometry/connectivity','volume':'independently regenerate and measure under same options; all generations count toward shared limit','insufficient_budget':'FAIL verification; never exceed shared cap'},
 'synthetic':{'scales':[.5,1.,4.],'shapes':['ellipse','mild_irregular_convex_polygon'],'normalized_quality_tolerance':1e-10,'synthetic_length_units':'arbitrary, scale and volume target multiplied together'},
 'contract_sha256':b['contract_sha256'],'baseline_mesh_sha256':b['mesh_sha256'],'baseline_recomputed_sha256':sha256(r/'baseline_recomputed.json')}
write_json(r/'acceptance_policy.json',policy)
lock={'timestamp':timestamp(),'policy_sha256':sha256(r/'acceptance_policy.json'),'contract_sha256':b['contract_sha256'],'baseline_mesh_sha256':b['mesh_sha256'],
 'source_surface_sha256':b['source_surface_sha256'],'baseline_recomputed_sha256':sha256(r/'baseline_recomputed.json'),'frozen_before_any_adaptive_search':True}
write_json(r/'freeze_lock.json',lock)
for name in ('acceptance_policy.json','freeze_lock.json'):shutil.copyfile(r/name,i/name)
print(json.dumps(lock,indent=2))
