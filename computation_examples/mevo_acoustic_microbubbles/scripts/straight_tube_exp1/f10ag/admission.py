"""Keep solver completion, model pairing and refinement scopes distinct."""
def time_windows(science):
    convergence = science['convergence']
    complete = [r for r in science['rows'] if r['kind'] == 'HF' and r['status'] == 'PASS']
    paired = {r['id']: r for r in science['paired'] if r['PASS']}
    rows = []
    for r in complete:
        if r['task'].get('variant'):
            continue
        checks = {q['id'].rsplit('_', 1)[-1]: q for q in convergence if q['parent'] == r['id']}
        refined = all(k in checks and checks[k]['PASS'] for k in ['DT05', 'FINE'])
        rows.append(dict(id=r['id'],state=r['state'],action=r['action'],N=r['N'],
                         completed_s=r['completed_s'],paired_pass=r['id'] in paired,
                         checks=checks,both_dt_and_cache_refined=refined,
                         scope='This exact initial state, action and completed window only'))
    return dict(rows=rows,
                max_new_HF_dt_and_cache_refined_s=max([r['completed_s'] for r in rows if r['both_dt_and_cache_refined']]+[0]),
                max_new_paired_and_refined_s=max([r['completed_s'] for r in rows if r['both_dt_and_cache_refined'] and r['paired_pass']]+[0]),
                historical_refinement='F10AF 50 ms X_PAIR convergence remains historical evidence, not a newly rerun AG refinement',
                limit_m=50e-9)
