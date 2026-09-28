"""Exact read-only extraction of the validated native log parser."""
import re, math
NUMBER=r'[-+]?(?:\d*\.?\d+(?:[eE][-+]?\d+)?|nan|inf)'

ROW=re.compile(r'^\s*(?P<eq>NS|SS)\s+(?P<step>\d+)-(?P<nl>\d+)(?P<saved>s?)\s+(?P<elapsed>'+NUMBER+r')\s+[\[!]\s*(?P<ndb>[-+]?\d+)\s+(?P<nr1>'+NUMBER+r')\s+(?P<nr0>'+NUMBER+r')\s+(?P<lr>'+NUMBER+r')[\]!]\s+(?P<open>[\[!])\s*(?P<iterations>\d+)\s+(?P<ldb>[-+]?\d+)\s+(?P<pct>\d+)[\]!]\s*(?P<message>.*)$',re.I)

def parse_base_log(text,dt):
    rows=[];pending=[];unparsed=[];petsc_reasons=[];monitors=[];current=[];active_reason=None
    for number,line in enumerate(text.splitlines(),1):
        monitor=re.match(r'^\s*(\d+)\s+KSP unpreconditioned resid norm\s+('+NUMBER+r')\s+true resid norm\s+('+NUMBER+r')\s+\|\|r\(i\)\|\|/\|\|b\|\|\s+('+NUMBER+r')',line,re.I)
        if monitor:
            values=[float(v) for v in monitor.groups()[1:]]
            item={'line':number,'iteration':int(monitor[1]),'residual_norm':values[0] if math.isfinite(values[0]) else None,
                  'true_residual_norm':values[1] if math.isfinite(values[1]) else None,
                  'true_relative_residual':values[2] if math.isfinite(values[2]) else None,
                  'finite':all(math.isfinite(v) for v in values)}
            current.append(item)
        if 'ill-conditioned LHS' in line:
            pending.append({'type':'ILL_CONDITIONED_LHS','line':number,'text':line.strip()})
        if 'DIVERGED_' in line or 'CONVERGED_' in line:
            reason=re.search(r'(?:DIVERGED|CONVERGED)_\w+',line)
            iterations=re.search(r'iterations\s+(\d+)',line)
            active_reason={'line':number,'text':line.strip(),'reason':reason[0] if reason else None,
                           'iterations':int(iterations[1]) if iterations else None,'diverged':'DIVERGED_' in line}
            petsc_reasons.append(active_reason)
        match=ROW.match(line)
        if not match:
            if re.match(r'^\s*(NS|SS)\s+\d+-',line):unparsed.append({'line':number,'text':line})
            continue
        g=match.groupdict();step=int(g['step']);nl=int(g['nl'])
        floats={k:float(g[k]) for k in ('elapsed','nr1','nr0','lr')}
        finite=all(math.isfinite(v) for v in floats.values())
        success=g['open']=='[' and 'not converged' not in g['message'] and finite
        warnings=pending[:]
        if not success:warnings.append({'type':'LINEAR_NONCONVERGENCE' if finite else 'NONFINITE_RESIDUAL','line':number,'text':g['message']})
        rows.append({'step':step,'time_s':step*dt,'nonlinear_iteration':nl,'linear_solve_index':len(rows)+1,
                     'log_line':number,'saved_state':bool(g['saved']),'reported_elapsed_s':floats['elapsed'] if math.isfinite(floats['elapsed']) else None,
                     'linear_iterations':int(g['iterations']),'linear_converged':success,
                     'linear_residual_relative':floats['lr'] if math.isfinite(floats['lr']) else None,
                     'linear_residual_absolute':None,'residual_finite':finite,
                     'residual_trustworthy':finite and not bool(pending),
                     'nonlinear_Ri_over_R1':floats['nr1'] if math.isfinite(floats['nr1']) else None,
                     'nonlinear_Ri_over_R0':floats['nr0'] if math.isfinite(floats['nr0']) else None,
                     'linear_reported_db':int(g['ldb']),'nonlinear_reported_db':int(g['ndb']),
                     'linear_time_percent':int(g['pct']),'ill_conditioned_warning':bool(pending),'warnings':warnings})
        if active_reason is not None or current:
            rows[-1]['petsc_reason']=active_reason
            rows[-1]['petsc_monitor']=current[:]
            rows[-1]['linear_residual_absolute']=current[-1]['residual_norm'] if current else None
            if active_reason and active_reason['diverged']:rows[-1]['linear_converged']=False
            if any(not m['finite'] for m in current):
                rows[-1]['linear_converged']=False;rows[-1]['residual_finite']=False;rows[-1]['residual_trustworthy']=False
            monitors.extend(current);current=[];active_reason=None
        pending=[]
    return {'linear_solves':rows,'unparsed_rows':unparsed,'unassigned_warnings':pending,
            'petsc_reasons':petsc_reasons,
            'unassigned_petsc_monitor':current,'unassigned_petsc_reason':active_reason,
            'failed_linear_solves':sum(not r['linear_converged'] for r in rows),
            'ill_conditioned_warnings':text.count('ill-conditioned LHS'),
            'nonlinear_failure_messages':[l for l in text.splitlines() if 'number of nonlinear iterations' in l and 'exceeded' in l],
            'first_failed_step':next((r['step'] for r in rows if not r['linear_converged']),None)}


def parse_solver_log(text,dt):
    # Keep failed-attempt evidence separately; final nonlinear history refers to accepted corrections.
    # Raw logs are immutable. Nothing is removed from saved raw logs or all_petsc_attempt_reasons.
    original=parse_base_log(text,dt)
    pattern=r'SV13Q_BEGIN [^\n]*logical=(\d+) attempt=0 [^\n]*\n.*?SV13Q_END [^\n]*\n'
    def retained(match):
        block=match.group(0)
        end=block.rsplit('SV13Q_END ',1)[-1]
        logical=match.group(1)
        if 'healthy=0 ' in end and re.search(r'SV13Q_BEGIN logical='+logical+r' attempt=1 ',text):return ''
        return block
    accepted_text=re.sub(pattern,retained,text,flags=re.S)
    result=parse_base_log(accepted_text,dt)
    result['all_petsc_attempt_reasons']=original['petsc_reasons']
    result['recovered_attempts']=text.count('recovery=STALE_ILU_RECOVERED')
    return result
