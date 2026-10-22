"""Solver log parsing and numerical invariance checks."""
import csv
import gzip
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from vascular_validation.provenance import sha256,write_json
from vascular_validation.validation import require

ROOT=Path(__file__).resolve().parents[3]
REPORT=ROOT/'reports/solver_acceptance'
NUMBER=r'[-+]?(?:\d*\.?\d+(?:[eE][-+]?\d+)?|nan|inf)'
ROW=re.compile(r'^\s*(?P<eq>NS|SS)\s+(?P<step>\d+)-(?P<nl>\d+)(?P<saved>s?)\s+(?P<elapsed>'+NUMBER+r')\s+[\[!]\s*(?P<ndb>[-+]?\d+)\s+(?P<nr1>'+NUMBER+r')\s+(?P<nr0>'+NUMBER+r')\s+(?P<lr>'+NUMBER+r')[\]!]\s+(?P<open>[\[!])\s*(?P<iterations>\d+)\s+(?P<ldb>[-+]?\d+)\s+(?P<pct>\d+)[\]!]\s*(?P<message>.*)$',re.I)

def load(name,stage='solver_acceptance'):
    return json.loads((ROOT/'reports'/stage/(name+'.json')).read_text())

def parse_boundary(path,Q):
    lines=Path(path).read_text().splitlines()
    roles=next(s.split()[2:] for s in lines if s.startswith('step'))
    result=[]
    for text in lines:
        fields=text.split()
        if not fields or not fields[0].isdigit():continue
        values=list(map(float,fields[1:]))
        require(all(math.isfinite(x) for x in values),'Non-finite boundary integral')
        flow=dict(zip(roles,values[1:]))
        incoming=-flow['INLET'];total=sum(flow[n] for n in ('OUTLET_01','OUTLET_02','OUTLET_03'))
        result.append({'step':int(fields[0]),'time_s':values[0], 'signed_outward_flows_m3_s':flow,
                       'Q_in_m3_s':incoming,'Q_out_total_m3_s':total,
                       'epsilon_Q':abs(incoming-Q)/Q,'epsilon_mass':abs(total-incoming)/Q,
                       'signed_closure_over_Q':(total-incoming)/Q})
    return result

def parse_solver_log(text,dt):
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

def linear_gate(run,require_petsc=True):
    h=run['history'];rows=h['linear_solves']
    require(run['exit_code']==0 and not run.get('monitor_stop'),'Solver execution failed or was stopped by a hard gate')
    require(bool(rows) and not h['unparsed_rows'],'No complete linear solve history')
    require(not h['failed_linear_solves'] and not h['ill_conditioned_warnings'],'Linear convergence failed')
    require(all(r['residual_finite'] for r in rows),'Non-finite linear residual')
    if require_petsc:
        require(len(h['petsc_reasons'])==len(rows),'Missing or unmatched PETSc convergence reasons')
        require(not any(r['diverged'] for r in h['petsc_reasons']),'PETSc divergence overrides process exit zero')
        require(all(r.get('petsc_monitor') and all(m['finite'] for m in r['petsc_monitor']) for r in rows),'Missing or nonfinite PETSc residual monitor')
    return True

def nonlinear_gate(history,max_iterations=12,tolerance=1e-10,min_iterations=2):
    require(not history['nonlinear_failure_messages'],'Nonlinear failure warning')
    rows=history['linear_solves'];require(bool(rows),'No nonlinear iterations')
    final={r['step']:r for r in rows}
    for row in final.values():
        require(row['nonlinear_iteration']>=min_iterations,'Incomplete nonlinear step')
        values=[row['nonlinear_Ri_over_R0'],row['nonlinear_Ri_over_R1']]
        require(all(v is not None and math.isfinite(v) for v in values),'Nonfinite nonlinear residual')
        require(min(values)<=tolerance,'Nonlinear iteration cap is not convergence')
    return True

def validate_petsc_xml(path):
    ls=ET.parse(path).find('.//LS');require(ls is not None and ls.get('type')=='GMRES','Expected supported GMRES')
    require(all(c.tag in {'Linear_algebra','Max_iterations'} for c in ls),'Ignored or unsupported PETSc XML option')
    la=ls.find('Linear_algebra');require(la is not None and la.get('type')=='petsc','PETSc backend required')
    require(all(c.tag=='Preconditioner' for c in la),'Ignored PETSc linear algebra XML option')
    require(la.findtext('Preconditioner') in {'petsc-jacobi','petsc-rcs'},'Unsupported XML preconditioner')
    return True

def measured_mass_gate(measurement):
    require(measurement.get('source')=='actual_vtu_surface_integration','Mass must come from measured solution, not target flow')
    q=measurement['Q_target_m3_s'];incoming=measurement['Q_in_m3_s'];out=measurement['outlet_flows_m3_s']
    require(set(out)=={'OUTLET_01','OUTLET_02','OUTLET_03'} and q>0,'Incomplete measured flows')
    require(all(math.isfinite(v) for v in [q,incoming,*out.values()]),'Nonfinite measured flows')
    error=abs(sum(out.values())-incoming)/q
    require(math.isclose(error,measurement['epsilon_mass'],rel_tol=1e-12,abs_tol=1e-16),'Reported mass error is inconsistent with measured flows')
    require(error<=1e-6 and abs(incoming-q)/q<=1e-6,'Actual solution fails mass gate')
    return True

def require_file_hash(path,expected):
    require(Path(path).is_file() and sha256(path)==expected,'Frozen input content changed: '+str(path))
    return True

def require_frozen_inputs():
    baseline=json.loads(gzip.decompress((REPORT/'history_baseline.json.gz').read_bytes()))
    checked=[]
    for path,item in baseline['files'].items():
        if path.startswith(('inputs/','outputs/mesh_and_flow/solver_mesh/','configs/')) and 'sha256' in item:
            require_file_hash(ROOT/path,item['sha256']);checked.append(path)
    require(bool(checked),'Missing frozen input inventory')
    return checked

def linear_settings(xml_path):
    tree=ET.parse(xml_path);ls=tree.find('.//LS');la=ls.find('Linear_algebra')
    return {'solver':ls.get('type'),'backend':la.get('type'),'preconditioner':la.findtext('Preconditioner'),
            'assembly':la.findtext('Assembly','none'),'settings':{c.tag:c.text for c in ls if c.tag!='Linear_algebra'}}

def compare_production_xml(before,after):
    a=ET.parse(before).getroot();b=ET.parse(after).getroot()
    def without_ls(element):
        return (element.tag,tuple(sorted(element.attrib.items())),(element.text or '').strip(),tuple(without_ls(c) for c in element if c.tag!='LS'))
    require(without_ls(a)==without_ls(b),'Production geometry, mesh paths, physics, BC, initial condition or time settings changed outside LS')
    return True

def write_csv(path,rows):
    require(bool(rows),'Empty parsed table')
    keys=list(rows[0]);
    with Path(path).open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=keys);w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in row.items()})
