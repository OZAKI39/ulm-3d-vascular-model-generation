"""Authorized geometry build and a single zero-step LAMMPS initialization."""
import argparse
import json
import os
import resource
import subprocess
import time
import numpy as np
import pyvista as pv
from .config import ROOT,dump
from .y_bifurcation_sources import OUT,RESULTS,YSources,digest
from .y_mural_clot_geometry import build
from .y_surface_coordinates import validate_and_map
from .interface_model import write_typed_data
from .calibrate_failure import prelude
from .fragment_analysis import read_dump


def load_geometry():
    z=np.load(OUT/'CURVED_MURAL_CLOT_GEOMETRY.npz');a=np.load(OUT/'CURVED_WALL_ANCHORS.npz');b=np.load(OUT/'BOND_GROUPS.npz')
    n=len(z['X']);N=n+len(a['X'])
    return dict(X=np.vstack([z['X'],a['X']]),clot_count=n,anchors=np.arange(N)>=n,
                bonds=b['bonds_zero_based'],bond_types=b['bond_type'],surface=np.r_[z['exposed_surface'],np.zeros(len(a['X']),bool)],
                base=np.r_[z['attachment_surface'],np.zeros(len(a['X']),bool)],
                areas=np.r_[z['areas_m2'],np.zeros(len(a['X']))],normals=np.vstack([z['normals'],np.zeros_like(a['X'])]),
                info=json.loads((OUT/'GEOMETRY_SUMMARY.json').read_text()))


def geometry(configuration=None):
    started=time.perf_counter();s=YSources();audit=s.audit();assert audit['status']=='PASS'
    requested=json.loads(open(configuration).read()) if configuration else None
    g,c,solid,top,attach=build(s,requested)
    validate_and_map(s,g,c,solid,top,attach)
    s.spot_check(g['X'][np.flatnonzero(g['surface'])[::max(1,int(g['surface'].sum()/4))]][:4])
    # Obtain historical hemisphere node count from its actual saved request/result.
    old=json.loads((ROOT/'verification/lammps_contact_cascade_diagnosis_and_minimal_fix_v1/CONTACT_MODIFIED/result.json').read_text())
    info=json.loads((OUT/'GEOMETRY_SUMMARY.json').read_text())
    info['old_hemisphere_nominal_volume_m3']=old['V_free0'];info['old_hemisphere_clot_nodes']=old['clot_nodes']
    dump(OUT/'GEOMETRY_SUMMARY.json',info)
    dump(OUT/'GEOMETRY_PERFORMANCE.json',dict(wall_seconds=time.perf_counter()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
         work='One actual mesh assembly and necessary geometry checks; no mechanics or CFD time steps',CLOT_PRESENT_CFD_SOLVED='NOT_RUN'))


def run0(*,out=OUT,geometry_loader=load_geometry):
    OUT=out
    assert not (OUT/'LAMMPS_RUN0_CHECK.json').exists(),'Completed run 0 is not automatically repeated'
    g=geometry_loader();c=json.loads((OUT/'PLACEMENT_CONFIGURATION.json').read_text())
    exe=ROOT/'tools/lammps-contact-onset-v1-build/lmp'
    prov=json.loads((ROOT/'verification/lammps_contact_transfer_and_fragment_size_v1/CONTACT_VARIANT_PROVENANCE.json').read_text())
    assert digest(exe)==prov['executable_sha256']
    write_typed_data(OUT/'Y_CLOT_LAMMPS.data',g,c)
    text=prelude(g,c,1).replace('official unmodified stable_30Sep2026','existing MODIFIED_1P0DX; geometry initialization only').replace('read_data body.data','read_data Y_CLOT_LAMMPS.data').replace('pair_style bpm/peri\n','pair_style bpm/peri contact_distance_factor 1.0\n')
    for name,mask in [('exposed_surface',g['surface']),('attachment_surface',g['base'])]:
        text+='group '+name+' id '+' '.join(map(str,np.flatnonzero(mask)+1))+'\n'
    text+='''group wall_anchor type 2
fix fixed anchors setforce 0 0 0
compute topology all property/local batom1 batom2 btype
dump topology all local 1 run0_bonds.local index c_topology[1] c_topology[2] c_topology[3]
dump points all custom 1 run0_points.lammpstrj id type x y z vx vy vz fx fy fz d_vfrac d_damage
dump_modify points sort id format float %.17g
thermo 1
thermo_style custom step time atoms bonds pe ke
thermo_modify lost error lost/bond error flush yes
# No integration fix and no applied load. Exactly zero time steps.
run 0
'''
    (OUT/'in.y_clot_geometry_check.lammps').write_text(text)
    cmd=[str(exe),'-in','in.y_clot_geometry_check.lammps','-log','LAMMPS_RUN0.log','-screen','LAMMPS_RUN0.screen']
    dump(OUT/'RUN0_COMMAND.json',dict(argv=cmd,cwd=OUT,executable_sha256=digest(exe),expected_time_steps=0,contact_model='MODIFIED_1P0DX'))
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    started=time.perf_counter()
    p=subprocess.run(cmd,cwd=OUT,env=env,capture_output=True,text=True,timeout=120)
    (OUT/'RUN0_LAUNCHER.log').write_text(p.stdout+p.stderr)
    log=(OUT/'LAMMPS_RUN0.log').read_text();assert p.returncode==0 and 'ERROR:' not in log
    frames=list(read_dump(OUT/'run0_points.lammpstrj'));assert len(frames)==1 and frames[0][0]==0
    _,cols,z=frames[0];z=z[np.argsort(z[:,0])]
    assert np.isfinite(z).all() and len(z)==len(g['X'])
    assert np.array_equal(z[:,0],np.arange(1,len(z)+1))
    pos=z[:,[cols.index(k) for k in ['x','y','z']]];assert np.max(np.abs(pos-g['X']))<1e-15
    assert np.array_equal(z[:,1],np.where(g['anchors'],2,1))
    assert np.max(np.abs(z[:,[cols.index(k) for k in ['vx','vy','vz']]]))==0
    assert np.max(np.abs(z[:,cols.index('d_damage')]))==0
    assert np.allclose(z[:,cols.index('d_vfrac')],c['dx_m']**3,rtol=1e-14,atol=0)
    bf=list(read_dump(OUT/'run0_bonds.local'));assert len(bf)==1 and bf[0][0]==0
    b=bf[0][2][:,1:4].astype(int);b[:,:2]=np.sort(b[:,:2],axis=1)-1;order=np.lexsort((b[:,1],b[:,0]));b=b[order]
    assert np.array_equal(b[:,:2],g['bonds']) and np.array_equal(b[:,2],g['bond_types'])
    report=dict(status='PASS',exit_code=p.returncode,step=0,time_s=0.,time_steps_advanced=0,atoms=len(z),clot_nodes=g['clot_count'],
                wall_anchors=int(g['anchors'].sum()),bonds=len(b),bond_type_counts={str(k):int((b[:,2]==k).sum()) for k in (1,2,3)},
                coordinates_unchanged=True,velocities_zero=True,damage_zero=True,nodal_volumes_match=True,topology_matches=True,
                fixed_anchor_group_loaded=True,exposed_group_loaded=True,attachment_group_loaded=True,wall_seconds=time.perf_counter()-started,
                executable_sha256=digest(exe),maximum_initial_force_N=float(np.max(np.linalg.norm(z[:,[cols.index(k) for k in ['fx','fy','fz']]],axis=1))),
                note='read_data + run 0 only; expected no-time-integration warning; no fragmentation result',CLOT_PRESENT_CFD_SOLVED='NOT_RUN')
    dump(OUT/'LAMMPS_RUN0_CHECK.json',report);print(report,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['geometry','run0']);parser.add_argument('--configuration');args=parser.parse_args()
    if args.stage=='geometry':geometry(args.configuration)
    else:run0()
