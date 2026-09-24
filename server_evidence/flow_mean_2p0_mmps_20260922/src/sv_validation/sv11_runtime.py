"""Owned PETSc configuration and monitored execution, exclusively in SV1.1."""
import os, re, signal, subprocess, time
import xml.etree.ElementTree as ET
from pathlib import Path
from .sv11 import ROOT, REPORT, load, parse_solver_log
from .provenance import write_json

# These are PETSc 3.19.6 options read by the pinned KSPSetFromOptions interface.
OPTIONS={
    '-skip_petscrc':'', '-ksp_type':'gmres', '-ksp_pc_side':'right',
    '-ksp_norm_type':'unpreconditioned', '-ksp_rtol':'1e-10', '-ksp_atol':'1e-24',
    '-ksp_max_it':'2000', '-ksp_gmres_restart':'100',
    '-ksp_diagonal_scale':'', '-ksp_diagonal_scale_fix':'',
    '-pc_type':'asm', '-pc_asm_overlap':'2', '-sub_ksp_type':'preonly',
    '-sub_pc_type':'ilu', '-sub_pc_factor_levels':'2',
    '-ksp_monitor_true_residual':'', '-ksp_converged_reason':'', '-ksp_view':'',
    '-options_view':'', '-options_left':''}

def petsc_ls(tree):
    ls=tree.find('.//LS');ls.clear();ls.set('type','GMRES')
    la=ET.SubElement(ls,'Linear_algebra',type='petsc')
    ET.SubElement(la,'Preconditioner').text='petsc-jacobi'
    ET.SubElement(ls,'Max_iterations').text='2000'
    # The XML parser does not deliver tolerance/restart to PETSc. Do not add
    # misleading XML fields; all effective values come from the owned options.
    return tree

def run_solver(case,name,dt,stop_step=None):
    manifest=load('petsc_build_manifest');assert manifest['status']=='PASS'
    value=' '.join(k+(' '+v if v else '') for k,v in OPTIONS.items())
    env=dict(os.environ,PATH='/usr/bin:/bin',LD_LIBRARY_PATH=manifest['runtime_library_path'],
             OMP_NUM_THREADS='1',PETSC_OPTIONS=value,LIBGL_ALWAYS_SOFTWARE='1')
    for key in ('PETSC_OPTIONS_YAML','PETSC_OPTIONS_TABLE','PETSC_OPTIONS_FILE'):env.pop(key,None)
    # In 3.19.6 PETSC_OPTIONS is read after rc files, so an environment-only
    # -skip_petscrc cannot establish isolation. Require actual absence instead.
    rc_paths=[Path.home()/'.petscrc',case/'.petscrc',case/'petscrc']
    assert not any(p.exists() for p in rc_paths),'Unexpected ambient PETSc options file'
    log=ROOT/'logs/sv1_1'/(name+'.log');resource=log.with_name(name+'_resources.txt')
    assert not log.exists(), 'Never overwrite a solver attempt'
    if stop_step is not None:(case/'STOP_SIM').write_text(str(stop_step)+'\n')
    command=['/usr/bin/mpiexec','--bind-to','none','-n','4',manifest['executable'],'solver.xml']
    start=time.monotonic();stop=None;carry='';offset=0
    tree=ET.parse(case/'solver.xml')
    max_nl=int(tree.findtext('.//Add_equation/Max_iterations'))
    tol_nl=float(tree.findtext('.//Add_equation/Tolerance'))
    with log.open('w') as out:
        process=subprocess.Popen(['/usr/bin/time','-v','-o',str(resource),*command],cwd=case,env=env,
                                 stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
        while process.poll() is None:
            time.sleep(.5)
            with log.open() as reader:reader.seek(offset);chunk=reader.read();offset=reader.tell()
            joined=carry+chunk
            if '\n' in joined:complete,carry=joined.rsplit('\n',1)
            else:complete='';carry=joined
            if stop is None:
                parsed=parse_solver_log(complete,dt)
                bad=next((r for r in parsed['linear_solves'] if not r['linear_converged']),None)
                error=re.search(r'DIVERGED_\w+|PETSC ERROR|PC failed due to[^\n]*|(?<![a-z])(?:nan|inf)(?![a-z])',complete,re.I)
                nonlinear=parsed['nonlinear_failure_messages']
                nonlinear=nonlinear or [r for r in parsed['linear_solves'] if r['nonlinear_iteration']>=max_nl and
                                         (r['nonlinear_Ri_over_R0'] is None or r['nonlinear_Ri_over_R1'] is None or
                                          min(r['nonlinear_Ri_over_R0'],r['nonlinear_Ri_over_R1'])>tol_nl)]
                if bad or error or nonlinear:
                    stop={'reason':error.group() if error else ('LINEAR_NONCONVERGENCE' if bad else 'NONLINEAR_CONVERGENCE_FAIL'),
                          'row':bad,'elapsed_s':time.monotonic()-start}
                    # Official stop file is checked at timestep boundaries. For a
                    # PETSc API error use SIGTERM to avoid invalid downstream work.
                    (case/'STOP_SIM').write_text('0\n')
                    write_json(REPORT/(name+'_stop.json'),stop)
                    if error and ('PETSC ERROR' in error.group() or 'nan' in error.group().lower() or 'inf' in error.group().lower()):
                        os.killpg(process.pid,signal.SIGTERM)
                if parsed['unparsed_rows']:
                    write_json(REPORT/(name+'_unparsed.json'),parsed['unparsed_rows'])
            elif time.monotonic()-start-stop['elapsed_s']>60:
                os.killpg(process.pid,signal.SIGTERM)
        code=process.wait()
    elapsed=time.monotonic()-start;text=log.read_text(errors='replace')
    rss=re.search(r'Maximum resident set size \(kbytes\): (\d+)',resource.read_text())
    result={'command':command,'cwd':str(case),'exit_code':code,'elapsed_s':elapsed,
            'peak_rss_kib':int(rss[1]) if rss else None,'memory_scope':'GNU time maximum process RSS; not MPI aggregate',
            'mpi_ranks':4,'omp_num_threads':1,'log':str(log.relative_to(ROOT)),
            'resource_log':str(resource.relative_to(ROOT)),'PETSC_OPTIONS':value,'monitor_stop':stop,
            'ambient_options_files_checked_absent':[str(p) for p in rc_paths],
            'history':parse_solver_log(text,dt)}
    write_json(REPORT/(name+'_execution.json'),result)
    return result
