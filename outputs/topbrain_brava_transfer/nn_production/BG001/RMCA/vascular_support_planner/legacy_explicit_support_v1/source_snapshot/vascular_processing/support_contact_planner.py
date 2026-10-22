"""Deterministic interval coverage along branch arc, with transparent penalties."""
import numpy as np
import trimesh
from .support_inputs import save,table

def influence(point,cfg):
    # Smaller radius and tighter bends shorten admissible unsupported reach.
    base=cfg['critical_max_span_mm']/2
    curvature=point['local_curvature'];radius=point['local_radius']
    reach=base*min(1.,max(.55,radius/.65))/(1+min(.5,curvature*radius))
    return max(cfg['sample_step_mm'],reach)

def coverage_sets(raw,cfg):
    critical=[p for p in raw if p['risk_level']=='CRITICAL'];sets=[]
    for p in raw:
        reach=influence(p,cfg)
        covered={i for i,q in enumerate(critical) if q['branch_id']==p['branch_id'] and abs(q['arc_position']-p['arc_position'])<=reach+1e-9}
        # An initial unsupported island must have a tip at its lowest print
        # layer: distant arc coverage alone cannot start an island in mid-air.
        for i in list(covered):
            q=critical[i]
            if 'LOCAL_MINIMUM' in q['detectors_triggered'] and (abs(q['arc_position']-p['arc_position'])>.51 or p['xyz'][2]>q['xyz'][2]+.1):covered.remove(i)
        sets.append(covered)
    return critical,sets

def prune(raw,cfg):
    critical,sets=coverage_sets(raw,cfg);remaining=set(range(len(critical)));selected=[]
    while remaining:
        candidates=[]
        for i,p in enumerate(raw):
            new=sets[i]&remaining
            if not new or not p['surface_hit_valid']:continue
            # No hidden weighted total score: coverage, keep-out, radius,
            # curvature, stable ID are explicitly ordered.
            rank=(-len(new),p['distance_to_bifurcation']<cfg['bifurcation_keepout_mm'],1/max(p['local_radius'],1e-6),p['local_curvature'],p['point_id'])
            if cfg.get('selection_priority')=='vascular':
                rank=(p['distance_to_bifurcation']<cfg['bifurcation_keepout_mm'],1/max(p['local_radius'],1e-6),p['local_curvature'],-p.get('removal_accessibility',0),-len(new),p['point_id'])
            candidates.append((rank,i))
        if not candidates:break
        _,i=min(candidates);selected.append(i);remaining-=sets[i]
    # Remove redundant contacts, first removing keep-out and thin contacts.
    for i in sorted(selected,key=lambda j:(raw[j]['distance_to_bifurcation']>=cfg['bifurcation_keepout_mm'],raw[j]['local_radius'],raw[j]['point_id'])):
        others=set().union(*(sets[j] for j in selected if j!=i))
        if sets[i]<=others:selected.remove(i)
    # Coverage-preserving one-for-one exchanges can move a tip away from a
    # fork after the greedy pass without increasing the number of contacts.
    def penalty(i):
        p=raw[i];return (p['distance_to_bifurcation']<cfg['bifurcation_keepout_mm'],1/max(p['local_radius'],1e-6),p['local_curvature'],p['point_id'])
    for i in list(selected):
        others=set().union(*(sets[j] for j in selected if j!=i));required=sets[i]-others
        choices=[j for j in range(len(raw)) if raw[j]['surface_hit_valid'] and required<=sets[j] and penalty(j)<penalty(i)]
        if choices:selected[selected.index(i)]=min(choices,key=penalty)
    contacts=[]
    for n,i in enumerate(sorted(selected,key=lambda j:(raw[j]['branch_id'],raw[j]['arc_position']))):
        p=raw[i];inside=p['distance_to_bifurcation']<cfg['bifurcation_keepout_mm']
        contacts.append(dict(**p,contact_id=f'C{n:03}',tip_diameter=.6,coverage_radius_arc_mm=influence(p,cfg),critical_regions_covered=[critical[k]['point_id'] for k in sorted(sets[i])],
            thin_vessel_penalty=1/p['local_radius'],keepout_penalty=int(inside),curvature_penalty=p['local_curvature'],status='KEEP_OUT_EXCEPTION_CRITICAL_COVERAGE' if inside else 'SELECTED',
            keepout_reason='Removing this contact would leave critical requirements uncovered' if inside else None))
    return contacts,coverage(raw,contacts,cfg)

def coverage(raw,contacts,cfg):
    critical,sets=coverage_sets(raw,cfg);idx={p['point_id']:i for i,p in enumerate(raw)}
    covered=set().union(*(sets[idx[p['point_id']]] for p in contacts if p['point_id'] in idx)) if contacts else set()
    gaps=[];spacings=[]
    for branch in sorted({p['branch_id'] for p in critical}):
        points=sorted(p['arc_position'] for p in contacts if p['branch_id']==branch)
        spacings.extend(np.diff(points).tolist())
        uncovered=sorted(q['arc_position'] for i,q in enumerate(critical) if i not in covered and q['branch_id']==branch)
        if uncovered:
            splits=np.split(uncovered,np.flatnonzero(np.diff(uncovered)>cfg['sample_step_mm']*1.6)+1)
            gaps.extend(float(a[-1]-a[0]+cfg['sample_step_mm']) for a in splits)
    return dict(percent=100*len(covered)/len(critical) if critical else 100.,critical_requirements=len(critical),covered_requirements=len(covered),
        uncovered_point_ids=[q['point_id'] for i,q in enumerate(critical) if i not in covered],max_critical_unsupported_span=max(gaps,default=0.),
        contact_arc_spacing_mm=dict(minimum=min(spacings) if spacings else None,median=float(np.median(spacings)) if spacings else None,maximum=max(spacings) if spacings else None),
        estimated_arc_influence=True,actual_slicer_coverage_required=True)

def masks(mesh,contacts,tip_diameter=.6):
    """Discrete bottom patches; everything outside these is a support blocker.

    Boolean difference makes enforcer and blocker mutually exclusive. The
    blocker encloses upper surface, bifurcations, and already covered spans.
    These are modifier regions, never printable supports or vascular edits.
    """
    patches=[]
    for p in contacts:
        xyz=np.array(p['xyz']);r=max(tip_diameter/2,.25)
        disk=trimesh.creation.cylinder(radius=r,height=.7,sections=24);disk.apply_translation(xyz+[0,0,.10]);patches.append(disk)
    if not patches:raise ValueError('NO_SUPPORT_CONTACTS')
    enforcer=trimesh.boolean.union(patches,engine='manifold')
    # A closed blocker with enforcer holes suppresses automatic wraparound.
    blocker=trimesh.boolean.difference([mesh,enforcer],engine='manifold')
    return enforcer,blocker

def write_contacts(out,contacts,cov,cfg):
    table(out/'tables/support_contacts.csv',contacts);save(out/'planned_coverage.json',cov)
    save(out/'contact_plan.json',dict(contacts=contacts,coverage=cov,keepout_mm=cfg['bifurcation_keepout_mm'],
        selection_order=(['critical coverage is mandatory','keepout penalty','inverse local radius','curvature','top accessibility','new critical arc requirements covered','stable point ID'] if cfg.get('selection_priority')=='vascular' else ['new critical arc requirements covered','keepout penalty','inverse local radius','curvature','stable point ID']),
        pruning='Greedy branch-arc set cover followed by deletion of every redundant selected contact'))
