#!/usr/bin/env python3
"""Execute the frozen native continuation and monitor every solve and saved VTU."""
import csv,json,math,os,re,signal,statistics,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv12 import REPORT,OUTPUT,CONFIG,LOG,load,validate_frozen,validate_xml,interval_metrics,steady_gate,consecutive_passes,extension_eligible,checkpoint_audit
from sv_validation.sv11 import ROW,parse_solver_log,nonlinear_gate,parse_boundary,write_csv
from sv_validation.provenance import write_json,sha256,now
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.validation import require,integral_agreement

case=OUTPUT/'vascular_flow';policy=json.loads((CONFIG/'policy.json').read_text());dt=policy['dt_s']
build=load('petsc_build_manifest','sv1_1');options=json.loads((CONFIG/'petsc_options.json').read_text())
q=load('petsc_short_qc','sv1_1')['states'][-1]
measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',q['Q_target_m3_s'],policy['Umean_m_s'])
qc={'states':[],'intervals':[],'last_step':10,'steady_reached':False,'consecutive_passing_intervals':0}
previous_u=None;previous_q=None

def add_state(path):
    global previous_u,previous_q
    u,p=measure.read(path);state=measure.measure(u,p);state.pop('outlet_fractions',None)
    step=int(path.stem.rsplit('_',1)[1]);state.update(step=step,time_s=step*dt,path=str(path.relative_to(ROOT)),
        sha256=sha256(path),size_bytes=path.stat().st_size,source='actual_vtu_surface_integration',
        origin='SV1.1 native restart seed' if step==10 else 'SV1.2 native solver')
    if previous_q:qc['intervals'].append(interval_metrics(measure,u,previous_u,state,previous_q,policy['save_interval_steps']))
    previous_u=u;previous_q=state;qc['states'].append(state);qc['last_step']=step
    qc['steady_reached']=steady_gate(qc['intervals']);qc['consecutive_passing_intervals']=consecutive_passes(qc['intervals'])
    write_json(REPORT/'saved_state_qc.json',qc)
    rows=[{'step':s['step'],'time_s':s['time_s'],'Qtarget':s['Q_target_m3_s'],'Qin':s['Q_in_m3_s'],
           **s['outlet_flows_m3_s'],'Qout_total':s['Q_out_total_m3_s'],'epsilon_Q':s['epsilon_Q'],'epsilon_mass':s['epsilon_mass'],
           'velocity_finite':s['velocity_finite'],'pressure_finite':s['pressure_finite'],'wall_velocity_max_m_s':s['wall_velocity_max_m_s'],
           'source':s['source'],'origin':s['origin']} for s in qc['states']]
    write_csv(OUTPUT/'qc/mass_history.csv',rows)
    if qc['intervals']:write_csv(OUTPUT/'qc/steady_history.csv',qc['intervals'])
    print(f"VTU step {step}: epsilon_mass={state['epsilon_mass']:.6g}; consecutive steady intervals={qc['consecutive_passing_intervals']}",flush=True)
    return state

def run_block(target,extension=False):
    validate_frozen();validate_xml(case/'solver.xml',extension)
    assert not (case/'STOP_SIM').exists(),'Short-run STOP_SIM must not enter production'
    name='extension' if extension else 'production';log=LOG/(name+'.log');res=LOG/(name+'_resources.txt')
    assert not log.exists(),'Never overwrite a production attempt'
    for path in (Path.home()/'.petscrc',case/'.petscrc',case/'petscrc'):assert not path.exists(),'Unexpected PETSc ambient options file'
    env=dict(os.environ,PATH='/usr/bin:/bin',LD_LIBRARY_PATH=build['runtime_library_path'],
             PETSC_OPTIONS=options['PETSC_OPTIONS'],OMP_NUM_THREADS=str(options['OMP_NUM_THREADS']),LIBGL_ALWAYS_SOFTWARE='1')
    for key in ('PETSC_OPTIONS_YAML','PETSC_OPTIONS_TABLE','PETSC_OPTIONS_FILE'):env.pop(key,None)
    cmd=['/usr/bin/mpiexec','--bind-to','none','-n',str(options['mpi_ranks']),build['executable'],'solver.xml']
    start=time.monotonic();offset=0;carry='';pending=[];rows=[];stop=None;last_step=400 if extension else 10
    state_seen={s['step'] for s in qc['states']};stable={};last_progress=0;stream= (OUTPUT/'qc'/(name+'_linear_history.jsonl')).open('w')
    with log.open('w') as out:
        proc=subprocess.Popen(['/usr/bin/time','-v','-o',str(res),*cmd],cwd=case,env=env,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
        write_json(REPORT/'running.json',{'status':'RUNNING','pid':proc.pid,'started_utc':now(),'target_step':target,'name':name,'last_completed_solve':None})
        def request_stop(reason,detail):
            nonlocal stop
            if stop:return
            stop={'reason':reason,'detail':detail,'last_completed_solve':rows[-1] if rows else None,
                  'elapsed_s':time.monotonic()-start,'log_byte_offset':offset,'utc':now()}
            (case/'STOP_SIM').write_text('0\n');write_json(REPORT/(name+'_failure.json'),stop)
            print('SAFETY STOP:',reason,flush=True)
            if 'PETSC ERROR' in str(detail) or reason=='FIELD_INVALID':os.killpg(proc.pid,signal.SIGTERM)
        while True:
            alive=proc.poll() is None
            with log.open() as inp:inp.seek(offset);chunk=inp.read();offset=inp.tell()
            joined=carry+chunk
            if '\n' in joined:complete,carry=joined.rsplit('\n',1)
            else:complete='';carry=joined
            for line in complete.splitlines():
                pending.append(line)
                error=re.search(r'DIVERGED_\w+|PETSC ERROR|PC failed due to[^\n]*|(?<![a-z])(?:nan|inf)(?![a-z])',line,re.I)
                if error:request_stop('LINEAR_SOLVER_FAILURE',{'line':line,'pending_log':pending[-12:]})
                if 'Resetting restart flag to false' in line or "don't match" in line:request_stop('FROZEN_INPUT_CHANGED',line)
                if ROW.match(line):
                    parsed=parse_solver_log('\n'.join(pending),dt);row=parsed['linear_solves'][-1]
                    row['linear_solve_index']=len(rows)+1;rows.append(row);stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush();pending=[]
                    last_step=row['step']
                    if len(rows)==1 and last_step!=(401 if extension else 11):request_stop('FROZEN_INPUT_CHANGED','Native restart did not resume at expected step')
                    if not row['linear_converged'] or not row['residual_finite'] or not row.get('petsc_reason') or row['petsc_reason']['diverged']:
                        request_stop('LINEAR_SOLVER_FAILURE',row)
                    if row['nonlinear_iteration']>=12 and (row['nonlinear_Ri_over_R0'] is None or row['nonlinear_Ri_over_R1'] is None or min(row['nonlinear_Ri_over_R0'],row['nonlinear_Ri_over_R1'])>1e-10):
                        request_stop('NONLINEAR_CONVERGENCE_FAIL',row)
                elif re.match(r'^\s*NS\s+\d+-',line):request_stop('LINEAR_SOLVER_FAILURE','Unparseable solver summary: '+line)
            for path in sorted((case/'4-procs').glob('result_*.vtu')):
                step=int(path.stem.rsplit('_',1)[1])
                if step in state_seen:continue
                stat=path.stat();key=(stat.st_size,stat.st_mtime_ns)
                previous=stable.get(step);stable[step]=(key,time.monotonic() if previous is None or previous[0]!=key else previous[1])
                if not alive or (previous and previous[0]==key and time.monotonic()-previous[1]>=3):
                    try:add_state(path);state_seen.add(step)
                    except Exception as exc:request_stop('FIELD_INVALID',str(exc));state_seen.add(step)
            if time.monotonic()-last_progress>=15:
                write_json(REPORT/'running.json',{'status':'STOP_REQUESTED' if stop else 'RUNNING','pid':proc.pid,'target_step':target,
                   'name':name,'last_completed_solve':{k:v for k,v in rows[-1].items() if k!='petsc_monitor'} if rows else None,
                   'linear_solves_completed':len(rows),'saved_states':len(qc['states']),'elapsed_s':time.monotonic()-start})
                last_progress=time.monotonic()
            if not alive:break
            if stop and time.monotonic()-start-stop['elapsed_s']>60:os.killpg(proc.pid,signal.SIGTERM)
            time.sleep(.5)
        exit_code=proc.wait()
    stream.close();elapsed=time.monotonic()-start
    # Reparse the complete original log independently of the live stream.
    h=parse_solver_log(log.read_text(errors='replace'),dt)
    require(len(h['linear_solves'])==len(rows) and not h['unparsed_rows'],'Live and complete log histories disagree')
    final_by_step={r['step']:r for r in rows}
    nonlinear=[]
    for step,last in final_by_step.items():
        first=next(r for r in rows if r['step']==step)
        good=last['nonlinear_iteration']>=2 and last['residual_finite'] and min(last['nonlinear_Ri_over_R0'],last['nonlinear_Ri_over_R1'])<=1e-10
        nonlinear.append({'step':step,'time_s':step*dt,'iterations':last['nonlinear_iteration'],
                          'initial_residual':first['petsc_monitor'][0]['residual_norm'] if first.get('petsc_monitor') else None,
                          'final_residual':last['petsc_monitor'][0]['residual_norm'] if last.get('petsc_monitor') else None,
                          'Ri_over_R0':last['nonlinear_Ri_over_R0'],'Ri_over_R1':last['nonlinear_Ri_over_R1'],'converged':good})
    resource=res.read_text();rss=re.search(r'Maximum resident set size \(kbytes\): (\d+)',resource)
    wall=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)',resource)
    def seconds(s):
        value=0
        for part in s.split(':'):value=60*value+float(part)
        return value
    result={'command':cmd,'exit_code':exit_code,'elapsed_monotonic_s':elapsed,'GNU_wall_s':seconds(wall[1]) if wall else None,
            'peak_single_process_RSS_KiB':int(rss[1]) if rss else None,'memory_scope':'GNU time maximum process RSS, not MPI aggregate',
            'PETSC_OPTIONS':env['PETSC_OPTIONS'],'mpi_ranks':4,'OMP_NUM_THREADS':1,'log':str(log.relative_to(ROOT)),
            'resource_log':str(res.relative_to(ROOT)),'linear_solves':rows,'nonlinear_history':nonlinear,
            'linear_failures':h['failed_linear_solves']+sum(x['diverged'] for x in h['petsc_reasons'] if x['line']>(rows[-1]['log_line'] if rows else 0)),
            'nonlinear_failures':sum(not r['converged'] for r in nonlinear),'ill_conditioned_warnings':h['ill_conditioned_warnings'],
            'last_step':last_step,'target_step':target,'extension':extension,'stop':stop}
    # Use summary flags and all PETSc reasons without double counting matched failures.
    result['linear_failures']=max(h['failed_linear_solves'],sum(x['diverged'] for x in h['petsc_reasons']))
    write_json(REPORT/(name+'_execution.json'),result)
    if rows:write_csv(OUTPUT/'qc'/(name+'_linear_history.csv'),[{k:v for k,v in r.items() if k not in ('petsc_monitor','warnings')} for r in rows])
    if nonlinear:write_csv(OUTPUT/'qc'/(name+'_nonlinear_history.csv'),nonlinear)
    write_json(REPORT/'running.json',{'status':'FINISHED','exit_code':exit_code,'last_step':last_step,'stop':stop})
    return result

add_state(case/'4-procs/result_010.vtu')
require(sha256(case/'4-procs/stFile_last.bin')==load('restart_decision')['source_checkpoint_sha256'],'Checkpoint changed before launch')
run=run_block(400);runs=[run];extensions=0
if not run['stop'] and not run['linear_failures'] and not run['nonlinear_failures'] and run['last_step']==400 and not qc['steady_reached']:
    eligible=extension_eligible(run,qc);write_json(REPORT/'extension_decision.json',{'eligible':eligible,'used':eligible,'policy':policy['extension_eligibility']})
    if eligible:
        import xml.etree.ElementTree as ET
        cp=checkpoint_audit(case/'4-procs/stFile_last.bin',400,dt);write_json(REPORT/'extension_checkpoint.json',cp)
        tree=ET.parse(case/'solver.xml');tree.find('.//Number_of_time_steps').text='800';tree.write(CONFIG/'sv_flow_extension.xml',encoding='utf-8',xml_declaration=True)
        import shutil
        shutil.copyfile(CONFIG/'sv_flow_extension.xml',case/'solver.xml');extensions=1;runs.append(run_block(800,True))
final_run=runs[-1];failure=next((r['stop']['reason'] for r in runs if r['stop']),None)
if not failure and any(r['linear_failures'] for r in runs):failure='LINEAR_SOLVER_FAILURE'
if not failure and any(r['nonlinear_failures'] for r in runs):failure='NONLINEAR_CONVERGENCE_FAIL'
expected=list(range(10,(800 if extensions else 400)+1,policy['save_interval_steps']))
if not failure and ([s['step'] for s in qc['states']]!=expected or final_run['exit_code']!=0):failure='FIELD_INVALID'
if not failure and not qc['steady_reached']:failure='STEADY_NOT_REACHED'
final=qc['states'][-1]
if not failure and (not final['velocity_finite'] or not final['pressure_finite'] or not final['wall_noslip_pass']):failure='FIELD_INVALID'
if not failure and (final['epsilon_Q']>1e-6 or final['epsilon_mass']>1e-6):failure='MASS_BALANCE_FAIL'
if not failure:
    boundary=parse_boundary(case/'4-procs/B_NS_Velocity_flux.txt',measure.Q)
    row=next(r for r in boundary if r['step']==final['step']);native={r:row['signed_outward_flows_m3_s'][r] for r in final['signed_outward_boundary_flows_m3_s']}
    final['native_flux_difference_over_Q']=integral_agreement(native,final['signed_outward_boundary_flows_m3_s'],measure.Q)
    final['outlet_fractions']={r:value/final['Q_out_total_m3_s'] for r,value in final['outlet_flows_m3_s'].items()}
    final['accepted']=True;write_json(REPORT/'accepted_solution.json',final)
result={'status':'FAIL' if failure else 'NUMERICAL_PASS_PENDING_RELOAD_TESTS_REVIEW','reason':failure,'started_from':'native restart',
        'initial_step':10,'last_step':final_run['last_step'],'extension_used':extensions,'physical_time_s':final_run['last_step']*dt,
        'runs':['production']+(['extension'] if extensions else []),'accepted_solution_available':not bool(failure)}
write_json(REPORT/'flow_execution.json',result);print(json.dumps(result,indent=2),flush=True)
raise SystemExit(1 if failure else 0)
