"""Generate auditable SI decks; LAMMPS owns bond breaking and contact."""
from pathlib import Path
import numpy as np
from .config import dump
from .schema import validate

def write_data(path, geometry, material):
    X=geometry['X'];bonds=geometry['bonds'];pad=4*material['delta'];lo=X.min(0)-pad;hi=X.max(0)+pad
    with Path(path).open('w') as f:
        f.write(f'Equal-volume cubic clot body; generated independently of PD2\n\n{len(X)} atoms\n{len(bonds)} bonds\n\n1 atom types\n1 bond types\n\n')
        for axis,k in zip('xyz',range(3)):f.write(f'{lo[k]:.16g} {hi[k]:.16g} {axis}lo {axis}hi\n')
        f.write(f'\nMasses\n\n1 {material["nodal_mass_kg"]:.16g}\n\nAtoms # bond\n\n')
        for i,x in enumerate(X,1):f.write(f'{i} 1 1 {x[0]:.16g} {x[1]:.16g} {x[2]:.16g}\n')
        f.write('\nBonds\n\n')
        for i,(a,b) in enumerate(bonds,1):f.write(f'{i} 1 {a+1} {b+1}\n')

def force_conversion():
    return '''# Official BPM peri accumulates N/m^3: convert to nodal N BEFORE external forces.
variable si_fx atom fx*d_vfrac
variable si_fy atom fy*d_vfrac
variable si_fz atom fz*d_vfrac
fix si_force all setforce v_si_fx v_si_fy v_si_fz
'''

def write_lammps_case(request, geometry, material, load, directory):
    r=validate(request);p=Path(directory)
    if p.exists():raise FileExistsError(f'Case directory already exists: {p}')
    p.mkdir(parents=True);dump(p/'request.json',r);dump(p/'LOAD_PROVIDER.json',load['info']);dump(p/'material.json',material)
    write_data(p/'HEMISPHERE_LAMMPS.data',geometry,material)
    np.savez_compressed(p/'surface_forces.npz',surface_node_id=np.flatnonzero(geometry['surface'])+1,force_N=load['force_N'][geometry['surface']],source=r['source'])
    with (p/'surface_forces.in').open('w') as f:
        f.write('# Replaceable provider force table, SI N; material node IDs.\n')
        for i,F in enumerate(load['force_N'],1):
            if np.any(F):f.write(f'set atom {i} d_loadx {F[0]:.16g} d_loady {F[1]:.16g} d_loadz {F[2]:.16g}\n')
    base=' '.join(map(str,np.flatnonzero(geometry['base'])+1));surface=' '.join(map(str,np.flatnonzero(geometry['surface'])+1))
    steps=round((r['ramp_time_s']+r['hold_time_s'])/r['dt_s']);every=r['dump_every'];assert steps%every==0
    text=f'''# EXPLORATORY_ONLY — same geometry/material/time/load budget for all patterns.
units si
dimension 3
boundary m m m
atom_style bond
atom_modify map array
newton on off
special_bonds lj 0 1 1 coul 1 1 1
lattice sc {r['dx_m']:.16g}
read_data HEMISPHERE_LAMMPS.data extra/bond/per/atom 160 extra/special/per/atom 180
neighbor {2*r['dx_m']:.16g} bin
neigh_modify every 1 delay 0 check yes
comm_modify cutoff {r['delta_m']*2.5:.16g}
fix props all property/atom d_vfrac d_damage d_loadx d_loady d_loadz ghost yes
set group all d_vfrac {material['nodal_volume_m3']:.16g}
set group all d_loadx 0 d_loady 0 d_loadz 0
pair_style bpm/peri
pair_coeff * * {material['contact_micromodulus_N_m6']:.16g} {r['delta_m']:.16g}
bond_style bpm/peri store/local breaks {every} time id1 id2 x y z x/ref y/ref z/ref
bond_coeff 1 lps {r['K_Pa']:.16g} {r['G_Pa']:.16g} {r['delta_m']:.16g} {r['s00']:.16g} {r['alpha']:.16g}
group anchored id {base}
group free subtract all anchored
group surface id {surface}
include surface_forces.in
velocity all set 0 0 0 units box
{force_conversion()}
variable amp equal ternary(time<{r['ramp_time_s']:.16g},0.5*(1-cos(PI*time/{r['ramp_time_s']:.16g})),1.0)
variable load_fx atom v_amp*d_loadx
variable load_fy atom v_amp*d_loady
variable load_fz atom v_amp*d_loadz
fix stream surface addforce v_load_fx v_load_fy v_load_fz
fix drag free viscous {material['nodal_mass_kg']*r['numerical_drag_rate_s_inv']:.16g}
fix wall anchored setforce 0 0 0
fix integrate free nve
timestep {r['dt_s']:.16g}
compute raw_pe all pe
variable U_J equal c_raw_pe*{material['nodal_volume_m3']:.16g}
thermo {every*5}
thermo_style custom step time atoms bonds v_amp ke v_U_J
thermo_modify lost error lost/bond error flush yes format float %.12g
timer full
dump trajectory all custom {every} trajectory.lammpstrj id type x y z vx vy vz fx fy fz d_damage
dump_modify trajectory sort id format float %.15g
dump broken all local {every} broken_bonds.local index f_breaks[1] f_breaks[2] f_breaks[3] f_breaks[4] f_breaks[5] f_breaks[6] f_breaks[7] f_breaks[8] f_breaks[9]
dump_modify broken format float %.15g
compute topology all property/local batom1 batom2 btype
dump topology all local {steps} surviving_bonds.local index c_topology[1] c_topology[2] c_topology[3]
run {steps}
write_restart final.restart
'''
    (p/'in.lammps').write_text(text)
    return dict(path=p,request=r,geometry=geometry,material=material,load=load,steps=steps)
