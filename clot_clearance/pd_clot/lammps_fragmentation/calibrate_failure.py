"""At most four scalar critical-stretch solves per declared strength mapping."""
from pathlib import Path
import json,math
import numpy as np
from .config import dump
from .coherent_failure_card import OUTPUT
from .interface_model import calibration_geometry,write_typed_data
from .lammps_writer import force_conversion
from .runner import run_lammps
from .fragment_analysis import read_dump

def prelude(g,r,every):
    K=r['K_Pa'];G=r['G_Pa'];delta=r['delta_m'];dx=r['dx_m'];V=dx**3
    return f'''# EXPLORATORY_ONLY; official unmodified stable_30Sep2026
units si
dimension 3
boundary m m m
atom_style bond
atom_modify map array
newton on off
special_bonds lj 0 1 1 coul 1 1 1
lattice sc {dx:.17g}
read_data body.data extra/bond/per/atom 160 extra/special/per/atom 180
neighbor {2*dx:.17g} bin
neigh_modify every 1 delay 0 check yes
comm_modify cutoff {3*delta:.17g}
fix props all property/atom d_vfrac d_damage d_area d_refz ghost yes
set group all d_vfrac {V:.17g}
set group all d_area 0
variable initialz atom z
set group all d_refz v_initialz
group clot type 1
group anchors type 2
pair_style bpm/peri
pair_coeff * * {18*K/(math.pi*delta**4):.17g} {delta:.17g}
bond_style bpm/peri store/local breaks {every} time id1 id2 x y z x/ref y/ref z/ref
bond_coeff 1 lps {K:.17g} {G:.17g} {delta:.17g} {r['s00_bulk']:.17g} 0
bond_coeff 2 lps {K:.17g} {G:.17g} {delta:.17g} {r['s00_interface']:.17g} 0
bond_coeff 3 lps {K:.17g} {G:.17g} {delta:.17g} 1.0 0
velocity all set 0 0 0 units box
{force_conversion()}
'''

def dumps(every,steps,extra=''):
    return f'''dump trajectory all custom {every} trajectory.lammpstrj id type x y z vx vy vz fx fy fz d_damage {extra}
dump_modify trajectory sort id format float %.15g
dump broken all local {every} broken_bonds.local index f_breaks[1] f_breaks[2] f_breaks[3] f_breaks[4] f_breaks[5] f_breaks[6] f_breaks[7] f_breaks[8] f_breaks[9]
dump_modify broken format float %.15g
compute topology all property/local batom1 batom2 btype
dump topology all local {steps} surviving_bonds.local index c_topology[1] c_topology[2] c_topology[3]
run {steps}
write_restart final.restart
'''

def solve_calibration(r,kind,critical,trial,sbulk=None):
    p=OUTPUT/'calibration'/f'{kind}_{trial:02d}';p.mkdir(exist_ok=False)
    q={**r,'s00_bulk':critical if kind=='bulk' else sbulk,'s00_interface':critical if kind=='interface' else critical/2,'mpi_ranks':1}
    g=calibration_geometry(q,kind);write_typed_data(p/'body.data',g,q);dump(p/'request.json',q)
    steps=4000;every=40;strainmax=max(.8,8*critical);duration=strainmax/.01;dt=duration/steps;area=(8*q['dx_m'])**2
    text=prelude(g,q,every)+'variable zero equal 0.0\n'
    if kind=='bulk':
        mid=4*q['dx_m'];upper=' '.join(str(i+1) for i in np.flatnonzero(g['X'][:,2]>mid))
        text+=f'''group measure id {upper}
variable disp atom 0.01*time*(d_refz-{mid:.17g})
variable velo atom 0.01*(d_refz-{mid:.17g})
fix prescribed all move variable v_zero v_zero v_disp v_zero v_zero v_velo
'''
    else:
        text+=f'''group measure type 1
variable disp equal 0.01*time*{q['dx_m']:.17g}
variable velo equal 0.01*{q['dx_m']:.17g}
fix prescribed clot move variable v_zero v_zero v_disp v_zero v_zero v_velo
fix fixed anchors setforce 0 0 0
'''
    text+=f'''compute resultant measure reduce sum fz
variable stress equal -c_resultant/{area:.17g}
variable strain equal time*0.01
fix curve all ave/time 1 1 1 v_strain v_stress file nominal_stress.dat
timestep {dt:.17g}
thermo 100
thermo_style custom step time atoms bonds v_strain v_stress
thermo_modify lost error lost/bond error flush yes
timer full
'''+dumps(every,steps)
    (p/'in.lammps').write_text(text);case=dict(path=p,request=q,geometry=g,steps=steps);perf=run_lammps(case)
    a=np.loadtxt(p/'nominal_stress.dat',comments='#');peak=int(np.argmax(a[:,2]));events=[a for _,_,a in read_dump(p/'broken_bonds.local') if len(a)]
    row=dict(kind=kind,trial=trial,s00=critical,macro_peak_stress_Pa=float(a[peak,2]),nominal_strain_at_peak=float(a[peak,1]),
        cross_section_area_m2=area,broken_event_count=sum(len(a) for a in events),wall_seconds=perf['wall_seconds'],path=str(p),
        method='Affine laterally constrained uniaxial extension; cut-plane force / reference area' if kind=='bulk' else 'Rigid clot block separated normally from anchor layer; reaction / reference area',
        caveat='Strength-equivalent numerical cell; not reproduction of pure-shear/peeling geometry or toughness',prescribed_motion=True,free_inertial_dofs=0)
    dump(p/'CALIBRATION_RESULT.json',row);return row

def calibrate(r,kind,sbulk=None,trial_offset=0):
    target=[4000.,8000.] if kind=='bulk' else [1700.,3500.]
    # Predeclared scalar brackets, not an energy-to-stretch conversion.
    lo=0.;hi=.6 if kind=='bulk' else sbulk*(1-1e-6);trials=[]
    for i in range(1,5):
        critical=(lo+hi)/2;row=solve_calibration(r,kind,critical,i+trial_offset,sbulk);trials.append(row)
        print('CALIBRATION',kind,i+trial_offset,critical,row['macro_peak_stress_Pa'],flush=True)
        good=target[0]<=row['macro_peak_stress_Pa']<=target[1] and row['broken_event_count']>0
        if good:break
        if row['macro_peak_stress_Pa']<target[0]:lo=critical
        else:hi=critical
    result=dict(status='PASS' if good else 'FAIL',target_strength_range_Pa=target,s00=critical,
        calibrated_macro_strength_Pa=row['macro_peak_stress_Pa'],trials=trials,solve_count=len(trials),max_solves=4,
        trial_offset=trial_offset,anchor_layers=r.get('anchor_layers',1),
        method='Bounded scalar bisection; freeze first in-range solution',scope='EXPLORATORY_ONLY')
    dump(OUTPUT/(kind.upper()+'_CALIBRATION.json'),result)
    if not good:raise RuntimeError(f'{kind} strength calibration failed within four solves')
    return result
