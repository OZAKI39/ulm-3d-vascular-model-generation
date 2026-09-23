"""Dimensionless, bounded port search and volume feedback. No PDE objects."""
import math
import numpy as np
from .cap_remesh import triangle_quality
from .mesh_qc import quantiles


def rim_distance(points,polygon):
    points=np.atleast_2d(points);edges=np.roll(polygon,-1,axis=0)-polygon
    delta=points[:,None,:]-polygon[None,:,:]
    t=np.clip(np.sum(delta*edges[None,:,:],axis=-1)/np.sum(edges*edges,axis=1),0,1)
    return np.linalg.norm(delta-t[...,None]*edges,axis=-1).min(axis=1)


def size_function(distance,h_rim,H,grading_slope):
    if h_rim<=0 or H<=0 or grading_slope<=0: raise ValueError('Positive size parameters required')
    return np.minimum(H,h_rim+grading_slope*np.asarray(distance))


def surface_quality(points,triangles,policy):
    try:q,ratio=triangle_quality(points,triangles)
    except ValueError as exc:return {'status':'FAIL','reason':str(exc),'triangle_count':len(triangles)}
    stats=quantiles(q);g=policy['surface_quality']
    checks={'finite_positive':bool(np.isfinite(q).all() and (q>0).all()),'low_count':bool(np.count_nonzero(q<g['low_threshold'])==0),
            'P5':stats['P5']>=g['P5_min'],'median':stats['median']>=g['median_min']}
    return {'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'q_tri':stats,'low_count':int(np.count_nonzero(q<g['low_threshold'])),'triangle_count':len(triangles)}


def choose_surface(trials):
    good=[t for t in trials if t['geometry_status']=='PASS' and t['quality']['status']=='PASS']
    return min(good,key=lambda t:(t['triangle_count'],-t['H'])) if good else None


def surface_search(evaluate,h_volume,h_rim,R_eq,policy):
    """Actual evaluations only, bounded growth/shrink followed by log bisection."""
    trials=[];cfg=policy['search'];limit=policy['limits']['maximum_surface_trials_per_port']
    factor=cfg['search_factor'];H=h_volume;direction=None;hpass=None;hfail=None;refinements=0
    reason='INITIAL_VOLUME_TARGET';termination='MAX_SURFACE_TRIALS_REACHED'
    while len(trials)<limit:
        record=evaluate(H,len(trials),reason);record.update(H=float(H),search_reason=reason)
        trials.append(record)
        if record['geometry_status']!='PASS':termination='SURFACE_GEOMETRY_FAILED';break
        passed=record['quality']['status']=='PASS'
        if passed:hpass=H if hpass is None else max(hpass,H)
        else:hfail=H if hfail is None else min(hfail,H)
        if direction is None:direction='grow' if passed else 'shrink'
        if hpass is not None and hfail is not None:
            if hpass>=hfail:termination='NONMONOTONE_QUALITY_BRACKET';break
            if refinements>=cfg['maximum_bracket_refinements']:termination='BRACKET_REFINED';break
            H=math.sqrt(hpass*hfail);refinements+=1;reason='LOG_BRACKET_BISECTION';continue
        if direction=='grow':
            if H>=cfg['maximum_H_over_R_eq']*R_eq:termination='GEOMETRY_SCALE_LIMIT';break
            H=min(H*factor,cfg['maximum_H_over_R_eq']*R_eq);reason='QUALITY_PASS_TRY_COARSER'
        else:
            if H<=h_rim:termination='NO_PASS_AT_RIM_SCALE';break
            H=max(H/factor,h_rim);reason='QUALITY_FAIL_TRY_FINER'
    chosen=choose_surface(trials)
    return {'status':'PASS' if chosen is not None and termination!='SURFACE_GEOMETRY_FAILED' else 'FAIL','trials':trials,
            'selected_trial':chosen['trial'] if chosen else None,'selected_H':chosen['H'] if chosen else None,'termination_reason':termination,
            'refinement_evaluations':sum(t['search_reason']=='LOG_BRACKET_BISECTION' for t in trials)}


def volume_acceptance(measured,baseline,policy,loadable):
    q=measured['quality'];b=baseline['quality'];p=measured['proxy'];bp=baseline['proxy'];g=policy['volume_quality'];cost=policy['cost']
    tol=g['quantile_roundoff_tolerance'];quality={'cap_low':q['cap_adjacent_below_0_1']<=math.floor(g['cap_low_fraction']*b['cap_adjacent_below_0_1']),
      'total_low':q['total_below_0_1']<=math.floor(g['total_low_fraction']*b['total_below_0_1']),
      **{k:q['min_sicn'][k]+tol>=b['min_sicn'][k] for k in ('P1','P5','median')}}
    ratios={'C_P2':p['N_P2_velocity_proxy']/bp['N_P2_velocity_proxy'],'C_tetra':p['N_tetra']/bp['N_tetra']}
    costs={k:v<=cost[k+'_max'] for k,v in ratios.items()}
    validity=all(v==0 for v in measured['validity'].values()) and measured['topology']['connected_fluid_components']==1 and loadable
    return {'status':'PASS' if validity and all(quality.values()) and all(costs.values()) else 'FAIL',
     'validity_pass':bool(validity),'quality_pass':all(quality.values()),'quality_checks':quality,'cost_pass':all(costs.values()),'cost_checks':costs,**ratios}


def feedback(iteration,measured,acceptance,policy):
    """Stop on feasibility/cost; otherwise modify exactly one largest contributor."""
    limit=min(policy['limits']['maximum_volume_iterations'],policy['limits']['maximum_total_volume_meshes'])
    def stop(decision,reason):return {'decision':decision,'termination_reason':reason,'changed_port':None,'stop':True}
    if iteration>=limit:return stop('STOP_MAX_ITER','MAX_ADAPTIVE_ITERATIONS_REACHED')
    if not acceptance['validity_pass']:return stop('STOP_INVALID_VOLUME','VOLUME_VALIDITY_FAILED')
    if not acceptance['cost_pass']:return stop('STOP_COST','FEM_COST_BUDGET_EXCEEDED')
    if acceptance['quality_pass']:return stop('ACCEPT','FIRST_FEASIBLE_ACCEPTED')
    if iteration+1>=limit:return stop('STOP_MAX_ITER','MAX_ADAPTIVE_ITERATIONS_REACHED')
    if acceptance['quality_checks']['cap_low']:return stop('STOP_NOT_CAP_CORRECTABLE','NOT_CAP_CORRECTABLE')
    counts=measured['quality']['low_quality_nearest_boundary_counts']
    ports=sorted((name,count) for name,count in counts.items() if name!='WALL')
    name,count=min(ports,key=lambda item:(-item[1],item[0]))
    if count<=0:return stop('STOP_NOT_CAP_CORRECTABLE','NOT_CAP_CORRECTABLE')
    return {'decision':'REFINE_'+name,'changed_port':name.lower(),'stop':False,'termination_reason':None,
            'reason':'Cap low-count gate failed; refine only the port with the largest low-count contribution; lexical tie break',
            'H_divisor':policy['search']['volume_feedback_refinement_factor']}


def select_iteration(records):
    eligible=[r for r in records if r['surface_status']=='PASS' and r['acceptance']['status']=='PASS']
    return min(eligible,key=lambda r:(r['volume']['proxy']['N_P2_velocity_proxy'],r['volume']['proxy']['N_tetra'],r['volume']['quality']['cap_adjacent_below_0_1']))['iteration'] if eligible else None
