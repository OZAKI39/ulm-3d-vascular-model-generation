"""Separate rejected failed attempts from accepted, independently checked corrections.

No tolerance relaxation: a retry is acceptable only when original RHS/zero guess
are verified and its fresh solve actually converges before any NS correction row.
"""
import re
from flow_solver_support.flow_parser import parse_solver_log

def classify(text,dt):
    parsed=parse_solver_log(text,dt)
    blocks=[]
    for match in re.finditer(r'SV13Q_BEGIN ([^\n]*)\n.*?SV13Q_END ([^\n]*)\n',text,re.S):
        start=dict(re.findall(r'(\w+)=(\S+)',match[1]));end=dict(re.findall(r'(\w+)=(\S+)',match[2]))
        blocks.append(dict(begin=start,end=end,start=match.start(),stop=match.end()))
    unsafe=[];recovered=[]
    for block in blocks:
        end=block['end']
        if end.get('healthy')=='1':continue
        logical=end.get('logical');attempt=end.get('attempt')
        retry=next((b for b in blocks if b['start']>block['stop'] and b['begin'].get('logical')==logical and b['begin'].get('attempt')=='1'),None)
        between=text[block['stop']:retry['start']] if retry else ''
        marker=re.search(r'SV13Q_RECOVERY logical='+re.escape(logical or '')+r'\b[^\n]*original_rhs_restored=1\b[^\n]*initial_guess_zero=1\b',between)
        good=bool(attempt=='0' and retry and marker and retry['end'].get('healthy')=='1' and retry['end'].get('recovery')=='STALE_ILU_RECOVERED' and not re.search(r'^\s*NS\s+\d+-',between,re.M))
        detail=dict(logical=logical,failed_attempt=end,original_RHS_and_zero_guess_verified=bool(marker),retry=retry['end'] if retry else None,no_NS_correction_before_retry=not bool(re.search(r'^\s*NS\s+\d+-',between,re.M)))
        (recovered if good else unsafe).append(detail)
    attempts=parsed['all_petsc_attempt_reasons'];failed_attempts=sum(bool(a['diverged']) for a in attempts)
    solves=parsed['linear_solves'];ratios=[]
    for solve in solves:
        monitors=solve.get('petsc_monitor',[])
        if monitors:ratios.append(monitors[-1]['true_residual_norm']/max(1e-24,1e-10*monitors[0]['true_residual_norm']))
    finite=all(s['residual_finite'] for s in solves)
    fatal=re.findall(r'NaN|nan\b|Segmentation fault|FRESH_RETRY_FAILED',text)
    accepted=bool(solves and not unsafe and failed_attempts==len(recovered) and len(recovered)==parsed['recovered_attempts'] and not parsed['unparsed_rows'] and not parsed['failed_linear_solves'] and not parsed['nonlinear_failure_messages'] and finite and not fatal and ratios and max(ratios)<1.01)
    return dict(accepted=accepted,accepted_linear_corrections=len(solves),all_KSP_attempts=len(attempts),failed_KSP_attempts=failed_attempts,verified_recovered_attempts=recovered,unrecovered_or_unsafe_attempts=unsafe,failed_accepted_linear_corrections=parsed['failed_linear_solves'],unparsed_rows=len(parsed['unparsed_rows']),nonlinear_failure_messages=parsed['nonlinear_failure_messages'],fatal_markers=fatal,max_accepted_true_residual_over_criterion=max(ratios) if ratios else None,criteria='Original registered tolerances and <1.01 monitor numerical-comparison allowance; every failed attempt must restore original RHS, use zero guess, and freshly converge before any NS correction is accepted.')
