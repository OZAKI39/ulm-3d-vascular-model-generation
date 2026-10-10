"""One gated contact-onset candidate and two frozen-strength compatibility checks."""
from pathlib import Path
import json,subprocess,os,hashlib,shutil,time
import numpy as np
from .config import ROOT,dump
from . import runner
from .lammps_writer import force_conversion
from .fragment_analysis import read_dump
from .interface_model import calibration_geometry
from .contact_counterfactual import O,OLD,run,intervene

EXE=ROOT/'tools/lammps-contact-onset-v1-build/lmp'

def smoke():
    """Four stationary particles, run 0 only: force sign, cutoff, MPI, restart."""
    root=O/'checks/contact_style_smoke';root.mkdir(exist_ok=False)
    q=json.loads((O/'FROZEN_BASELINE.json').read_text())['request'];dx=q['dx_m'];delta=q['delta_m'];vol=dx**3;k=18*q['K_Pa']/(np.pi*delta**4)
    prefix=f'''units si
dimension 3
boundary f f f
atom_style bond
newton on off
processors 2 1 1
lattice sc {dx:.17g}
region box block -3 4.8 -2 2 -2 2 units lattice
create_box 1 box
create_atoms 1 single 0 0 0 units lattice
create_atoms 1 single 0.9 0 0 units lattice
create_atoms 1 single 3 0 0 units lattice
create_atoms 1 single 4.1 0 0 units lattice
mass 1 {q['rho_kg_m3']*vol:.17g}
fix props all property/atom d_vfrac ghost yes
set group all d_vfrac {vol:.17g}
'''
    tail=f'''neighbor {2*dx:.17g} bin
comm_modify cutoff {3*delta:.17g}
{force_conversion()}
thermo 1
dump snapshot all custom 1 forces.dump id x y z fx fy fz
dump_modify snapshot sort id format float %.17g
run 0
write_restart smoke.restart
'''
    specifications=[('official_default',runner.LMP,'',1.35),('variant_default',EXE,'',1.35),('variant_1p0',EXE,' contact_distance_factor 1.0',1.0)]
    outputs={};checks=[]
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    def execute(name,executable,text):
        p=root/name;p.mkdir();(p/'in.lammps').write_text(text)
        cmd=['mpirun','--bind-to','core','-np','2',str(executable),'-in','in.lammps','-log','log.lammps','-screen','screen.lammps']
        with (p/'launcher.log').open('w') as f:subprocess.run(cmd,cwd=p,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        a=list(read_dump(p/'forces.dump'))[-1][2];return a[np.argsort(a[:,0])][:,4:7]
    for name,exe,option,factor in specifications:
        F=execute(name,exe,prefix+'pair_style bpm/peri'+option+f'\npair_coeff * * {k:.17g} {delta:.17g}\n'+tail)
        expected=np.zeros((4,3))
        for i,j,r in [(0,1,.9*dx),(2,3,1.1*dx)]:
            f=15*k*vol**2/delta*max(factor*dx-r,0);expected[i,0]=-f;expected[j,0]=f
        assert np.allclose(F,expected,rtol=1e-12,atol=1e-18),(name,F,expected)
        assert np.max(abs(F.sum(0)))<1e-18
        outputs[name]=F;checks.append(dict(name=name,run_steps=0,atoms=4,MPI_ranks=2,force_matches_analytic=True,net_force_zero=True,max_force_error_N=float(abs(F-expected).max())))
    assert np.array_equal(outputs['official_default'],outputs['variant_default'])
    for previous in ['official_default','variant_1p0']:
        header=f'''units si
atom_style bond
newton on off
read_restart ../{previous}/smoke.restart
lattice sc {dx:.17g}
fix props all property/atom d_vfrac ghost yes
set group all d_vfrac {vol:.17g}
'''
        F=execute('restart_'+previous,EXE,header+tail);assert np.array_equal(F,outputs[previous])
        checks.append(dict(name='restart_'+previous,run_steps=0,atoms=4,MPI_ranks=2,forces_equal_original=True))
    result=dict(status='PASS',scope='Static code-level force smoke only; no integration, collision or benchmark',default_native_forces_bitwise_equal=True,nondefault_restart_factor_preserved=True,old_default_restart_readable=True,checks=checks)
    dump(O/'checks/CONTACT_STYLE_SMOKE.json',result);print('CONTACT_STYLE_SMOKE PASS',flush=True)

def compatibility(kind):
    source=OLD/'calibration'/('bulk_02' if kind=='bulk' else 'interface_08');p=O/'calibration_compatibility'/kind;p.mkdir(parents=True,exist_ok=False)
    shutil.copy2(source/'body.data',p/'body.data');r=json.loads((source/'request.json').read_text());r.update(contact_model='MODIFIED_BPM_PERI',compatibility_only=True);dump(p/'request.json',r)
    text=intervene((source/'in.lammps').read_text(),'MODIFIED_BPM_PERI');(p/'in.lammps').write_text(text);g=calibration_geometry(r,kind)
    saved=runner.LMP
    try:runner.LMP=EXE;perf=runner.run_lammps(dict(path=p,request=r,geometry=g,steps=4000))
    finally:runner.LMP=saved
    perf['scope']='One frozen-strength calibration compatibility check';dump(p/'PERFORMANCE.json',perf)
    curve=np.loadtxt(p/'nominal_stress.dat');i=int(np.argmax(curve[:,2]));limits=[4000,8000] if kind=='bulk' else [1700,3500];old=json.loads((source/'CALIBRATION_RESULT.json').read_text())
    result=dict(kind=kind,status='PASS' if limits[0]<=curve[i,2]<=limits[1] else 'FAIL',old_peak_Pa=old['macro_peak_stress_Pa'],new_peak_Pa=float(curve[i,2]),strain_at_peak=float(curve[i,1]),target_strength_range_Pa=limits,bulk_s00=r['s00_bulk'],interface_s00=r['s00_interface'],input_source=str(source),input_change='pair_style bpm/peri contact_distance_factor 1.0 only',recalibrated=False,solves=1,runtime=perf)
    dump(p/'COMPATIBILITY.json',result);(O/(kind.upper()+'_STRENGTH_COMPATIBILITY.md')).write_text(f"# {kind} strength compatibility\n\nStatus: {result['status']}. Original positive peak {old['macro_peak_stress_Pa']} Pa; candidate peak {result['new_peak_Pa']} Pa; frozen target {limits} Pa. Exactly one copied prescribed calibration input, only contact option changed. All s00 and other cell/load settings unchanged; no refit. Signed curve retained in calibration_compatibility/{kind}/nominal_stress.dat.\n")
    print('COMPATIBILITY',kind,result['status'],result['new_peak_Pa'],flush=True);return result

def execute_candidate():
    decision=json.loads((O/'STAGE_A_FROZEN_RESULT.json').read_text());assert decision['stage_B_allowed'] and decision['contact_role']=='PRIMARY_CONTRIBUTOR_SUPPORTED'
    smoke();run('MODIFIED_BPM_PERI',EXE)
    checks={kind:compatibility(kind) for kind in ['bulk','interface']};dump(O/'CALIBRATION_COMPATIBILITY.json',checks)
