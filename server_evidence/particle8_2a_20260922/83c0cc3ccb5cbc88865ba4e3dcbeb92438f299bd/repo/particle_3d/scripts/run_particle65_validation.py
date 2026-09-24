#!/usr/bin/env python3
"""Execute real physics and preserve raw data. Run before rendering the report."""
from pathlib import Path
import sys,json,tempfile
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.particle65_cases import *
from particle_3d.particle65_validation import static_scans,bridge_case
from particle_3d.particle3_cases import write_json,write_rows
from particle_3d.nearfield_regularization import contract
from particle_3d.nearfield_handoff import require_admissible_initial,InitialBelowContinuumHandoff
from particle_3d.particle6_checkpoint import frozen_provenance
REPORT=PACKAGE/'reports/particle6_5';DATA=REPORT/'data'
REAL_DT=4*1.0016153033502137e-5;REAL_HORIZON=500*REAL_DT


def save(name,value,rows=None):
    write_json(DATA/(name+'.json'),value)
    if rows is not None:write_rows(DATA/(name+'.csv'),rows)
    elif isinstance(value,list):write_rows(DATA/(name+'.csv'),value)
    else:write_rows(DATA/(name+'.csv'),[value])
    print('saved',name,flush=True)


def trajectory(name,case,dt,horizon,**kw):
    cached=DATA/(name+'.json')
    if '--reuse-trajectories' in sys.argv and cached.exists():
        d=json.loads(cached.read_text());assert d['summary']['dt_s']==dt and d['summary']['horizon_s']==horizon
    else:d=run_case(case,dt,horizon,**kw)
    if kw.get('old'):
        for r in d['interactions']:r['interaction_state_role']='V1_BOUNDARY_DIAGNOSTIC_ON_HISTORICAL_P5_STATE'
        d['summary']['parameter_role']='HISTORICAL_P5_DEFAULT; V1_SCALES_FOR_COMPARISON_ONLY'
    save(name,d,d['interactions'])
    write_rows(DATA/(name+'_positions.csv'),d['positions']);write_rows(DATA/(name+'_intervals.csv'),d['ledger'])
    print(d['summary'],flush=True);return d


def main():
    DATA.mkdir(parents=True,exist_ok=True)
    p=NearFieldRegularizationV1();frozen=frozen_provenance(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    historical={str(f.relative_to(REPO)):sha256(f) for f in sorted((PACKAGE/'reports/particle5').rglob('*')) if f.is_file()}
    p6=json.loads((PACKAGE/'reports/particle6/PARTICLE6_VALIDATION.json').read_text())
    assert all(sha256(REPO/k)==v for k,v in p6['source_sha256'].items())
    save('00_scope',dict(contract=contract(),frozen=frozen,p6_source_sha256=p6['source_sha256'],historical_p5_sha256=historical,
        viscosity_pa_s=MU,real_dt_s=REAL_DT,real_horizon_s=REAL_HORIZON,
        replay_window_choice='VALIDATION_ONLY; 20.032 ms includes return to WALL; same window and dt set for every floor',
        development_baseline='ALL_489_PASSED_BEFORE_V1_CODE'))
    chis=np.unique(np.r_[np.linspace(0,.06,601),.005,.01,.02,.03,.04,.05])
    save('01_activation',[dict(chi=x,w=p.activation_weight(x),derivative=p.activation_derivative(x)) for x in chis])
    save('02_lower',[dict(a_ref_m=a,**p.lower_handoff_gap(a)) for a in np.unique(np.r_[np.geomspace(.1e-6,10e-6,201),2e-6,.588015722765549e-6,1.2658269815352818e-6])])
    walls,pairs=static_scans();save('03_wall_static',walls);save('04_pair_static',pairs)
    summaries=[]
    for case,index in [('wall','05'),('pair','06')]:
        for div in [1,2,4]:
            r=trajectory(f'{index}_{case}_dt{div}',case,.004/div,.032);summaries.append(r['summary'])
        trajectory(f'{index}_{case}_old_p5',case,.004,.032,old=True)
    # Explicitly reject the original already sub-floor initial state without edits.
    shape,provider,wall,classifier,prov=real_fixture()
    original=json.loads((PACKAGE/'reports/particle5/data/10_initialization.json').read_text())['original_p4_initialization']['initial_state']['particles']
    old_shapes={r['particle_id']:Sphere(r['center_m'],r['radius_m']) for r in original}
    try:require_admissible_initial(old_shapes,wall,p,MU)
    except InitialBelowContinuumHandoff as error:rejected=error.record
    else:raise AssertionError('Invalid historical start should be rejected without movement')
    histpath=PACKAGE/'reports/particle5/data/10_real_two_mb_states.json';hist=json.loads(histpath.read_text());histrows=[]
    for s in hist:
        histrows.extend(dict(time_s=s['time_s'],particle_id=q['particle_id'],h_geom_m=q['wall_gap_m']) for q in s['particles'])
    save('07_historical_p5',dict(source=str(histpath.relative_to(REPO)),sha256=sha256(histpath),rows=histrows,
        minimum_h_geom_m=min(r['h_geom_m'] for r in histrows),original_two_mb_initial_rejection=rejected,
        replay_initialization=prov,comparison_scope='Historical two-MB context AND separate matched single-MB P5/V1 replay; do not conflate trajectories'),histrows)
    results={}
    for div in [1,2,4]:
        r=trajectory(f'07_real_v1_dt{div}','real',REAL_DT/div,REAL_HORIZON);summaries.append(r['summary']);results[(2.,div)]=r['summary']
    trajectory('07_real_matched_old_p5','real',REAL_DT,REAL_HORIZON,old=True)
    for nm in [1.5,3.]:
        for div in [1,2,4]:
            # Integer nanometre values are expressed explicitly, avoiding role validation by rounded arithmetic.
            floor=1.5e-9 if nm==1.5 else 3e-9
            r=trajectory(f'08_real_floor{nm:g}_dt{div}','real',REAL_DT/div,REAL_HORIZON,floor=floor);results[(nm,div)]=r['summary']
    sensitivity=[]
    for (nm,div),r in sorted(results.items()):
        base=results[(2.,div)];d=dict(r,dt_divisor=div,parameter_role='SENSITIVITY_ONLY',official_floor_m=2e-9,scientific_review='PENDING_USER_REVIEW')
        for key in ['handoff_time_s','minimum_h_geom_m','tangential_displacement_m','final_position_m','final_velocity_m_s','covered_time_s']:
            if r[key] is None or base[key] is None:d[key+'_absolute_difference_from_2nm']=None;d[key+'_relative_difference_from_2nm']=None;continue
            absolute=float(np.linalg.norm(np.asarray(r[key])-base[key]));denom=float(np.linalg.norm(base[key]))
            d[key+'_absolute_difference_from_2nm']=absolute;d[key+'_relative_difference_from_2nm']=absolute/denom if denom else (0. if absolute==0 else None)
        sensitivity.append(d)
    save('08_sensitivity',dict(rows=sensitivity,floors_m=[1.5e-9,2e-9,3e-9],same_real_state_fem_radius_wall_dt_set_horizon=True,
        metric_reference='same dt and 2 nm; vector differences use Euclidean norm',
        ordering={str(div):{key:dict(observed_ascending=[nm for nm in sorted([n for n in [1.5,2.,3.] if results[(n,div)][key] is not None],key=lambda n:results[(n,div)][key])],not_observed=[n for n in [1.5,2.,3.] if results[(n,div)][key] is None]) for key in ['minimum_h_geom_m','handoff_time_s','tangential_displacement_m']} for div in [1,2,4]},
        all_finite=all(r['finite'] for r in sensitivity),all_continuously_safe=all(r['all_accepted_above_lower'] for r in sensitivity),
        no_arbitrary_scientific_pass_threshold=True,scientific_sensitivity_review='PENDING_USER_REVIEW'),sensitivity)
    trajectory('09_large_step_refinement','wall',.032,.032)
    with tempfile.TemporaryDirectory(prefix='particle65_bridge_') as temp:
        b=bridge_case(Path(temp)/'checkpoint')
        # Preserve actual binary restart and verified metadata, replace only this stage's generated files.
        import shutil
        target=REPORT/'checkpoints/v1';target.mkdir(parents=True,exist_ok=True)
        for f in (Path(temp)/'checkpoint').iterdir():shutil.copyfile(f,target/f.name)
    save('10_bridge',b,[dict(step=r['step'],**{k:r[k] for k in ['R','b','U','J','standalone_eligible_pairs','lammps_candidates','filtered_nearfield_pairs','constraints_equal']}) for r in b['matrix_rows']])
    save('11_timestep',summaries)

if __name__=='__main__':main()
