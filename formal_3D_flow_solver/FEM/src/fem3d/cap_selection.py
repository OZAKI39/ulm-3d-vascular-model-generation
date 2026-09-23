"""Frozen Stage 1.5 candidate ordering; geometry rejection always takes priority."""
def require_volume_permission(surface):
    if surface['status']!='PASS' or not surface['volume_meshing_permitted']:
        raise ValueError('Geometry gate failed: volume meshing and mesh round-trip are forbidden')


def select_candidate(candidates,policy):
    eligible=[row for row in candidates if row['geometry_status']=='PASS' and row.get('volume_qc') is not None]
    def key(row):
        q=row['volume_qc']
        return (q['cap_adjacent_lt_0_1'],q['total_lt_0_1'],-q['P1'],-q['P5'],q['tetrahedra'],row['candidate'])
    ranked=sorted(eligible,key=key)
    if not ranked:
        return {'status':'FAIL','selected_candidate':None,'ranked_candidates':[],
                'reason':'No geometry-admissible candidate with a valid volume result; no winner, no replacement mesh'}
    winner=ranked[0];q=winner['volume_qc'];g=policy['quality']
    checks={'cap_tail':q['cap_adjacent_lt_0_1']<=g['cap_adjacent_lt_0_1_max'],
            'total_tail':q['total_lt_0_1']<=g['total_lt_0_1_max'],
            'P1':q['P1']>=g['P1_min'],'P5':q['P5']>=g['P5_min'],
            'median':q['median']>=g['median_min'],'validity':q['invalid_tetra']==0}
    return {'status':'CONDITIONAL PASS' if all(checks.values()) else 'FAIL',
            'selected_candidate':winner['candidate'],'quality_checks':checks,
            'ranked_candidates':[r['candidate'] for r in ranked],
            'ordering_values':{r['candidate']:list(key(r)) for r in ranked},
            'reason':'Lexicographic frozen priority: cap tail, total tail, P1, P5, tetra count; name breaks exact ties only'}
