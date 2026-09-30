"""File-import trajectory equivalence and small timestep refinement evidence."""
import copy,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pd_clot.runner import simulate
from pd_clot.output import write_json


def main():
    root=Path(__file__).resolve().parents[1];out=root/'verification'
    a=json.loads((root/'configs/analytic_static.json').read_text())
    b=json.loads((root/'configs/file_traction.json').read_text())
    ra=simulate(a,out/'analytic_static',root/'configs',quiet=True)
    rb=simulate(b,out/'file_traction',root/'configs',quiet=True)
    error=float(np.max(np.abs(ra['x']-rb['x'])))
    c=copy.deepcopy(b);c['streaming']['path']='../inputs/straight_pipe_poiseuille.vtu'
    c['name']='file_poiseuille_velocity_pressure';rc=simulate(c,out/'file_poiseuille',root/'configs',quiet=True)
    coarse=json.loads((root/'configs/straight_pipe.json').read_text())
    coarse['clot']['cells']=[6,4,4];coarse['clot']['origin_m']=[-.000375,-.00025,-.0011]
    coarse['simulation']['number_of_macro_steps']=2
    coarse['streaming']['bubble_center_m']=[0,0,-.0005]
    fine=copy.deepcopy(coarse);fine['simulation']['dt_s']/=2
    r0=simulate(coarse,out/'dt_coarse',quiet=True);r1=simulate(fine,out/'dt_half',quiet=True)
    displacement0=r0['x']-r0['cloud'].X;displacement1=r1['x']-r1['cloud'].X
    relative=float(np.linalg.norm(displacement1-displacement0)/np.linalg.norm(displacement1))
    damage_error=float(np.max(np.abs(r1['integrity']-r0['integrity'])))
    result=dict(status='PASS' if error<1e-12 and relative<.01 and damage_error<.005 else 'FAIL',
        static_analytic_vs_csv_max_position_error_m=error,velocity_pressure_import=rc['summary']['status'],
        dt_halving_displacement_relative_L2=relative,dt_halving_max_bond_integrity_error=damage_error,
        scope='One small 96-particle two-macro-step control; not mesh/horizon convergence or proof of local physics accuracy')
    write_json(out/'IMPORT_AND_DT_CHECK.json',result);print(json.dumps(result,indent=2))
    if result['status']!='PASS':raise SystemExit(1)


if __name__=='__main__':main()
