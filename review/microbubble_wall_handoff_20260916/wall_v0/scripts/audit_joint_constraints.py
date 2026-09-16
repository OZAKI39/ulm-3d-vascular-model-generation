from pathlib import Path
import json,subprocess,sys,numpy as np
R=Path(__file__).resolve().parents[1];W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['local_work']);sys.path.insert(0,str(R/'src'));from reference_wall_runtime import RuntimeReference
plan={'purpose':'Independent simultaneous wall + sphere-pair active-set static response','radius_m':9.683592065545495e-7,'wall_gap_over_a':.0001,'pair_gap_over_effective_radius':.0001,'dt_s':.001,'drive_left':[.005,0,-.005],'drive_right':[-.005,0,-.005],'does_not_launch_LAMMPS_or_advance_time':True};(R/'contracts/JOINT_CONSTRAINT_STATIC_PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
subprocess.run(['g++','-O2','-std=c++17','-I'+str(R/'src'),str(R/'src/audit_joint_constraints.cpp'),'-L'+str(W),'-lwall_math','-Wl,-rpath,'+str(W),'-o',str(W/'audit_joint_constraints')],check=True)
p=subprocess.run([str(W/'audit_joint_constraints'),str(R/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5')],capture_output=True,text=True);(R/'raw/JOINT_CONSTRAINT_STATIC_FIRST.json').write_text(p.stdout);assert p.returncode==0,p.stderr
out=json.loads(p.stdout);a=plan['radius_m'];sep=2*a+.0001*a/2;x=np.array([[-sep/2,0,1.0001*a],[sep/2,0,1.0001*a]])
c=json.loads((R/'cases/COMBINED_PAIR_WALL/CASE_CONTRACT.json').read_text());c['config'].update(field='WALL_ZERO_FIELD.h5',drive_x=0,drive_y=0,drive_z=0);ref=RuntimeReference(R,c)
original=ref.system
def system(x):
 A,b,active=original(x);b[0]+=.005*6*np.pi*.001*a;b[6]-=.005*6*np.pi*.001*a;b[2]-=.005*6*np.pi*.001*a;b[8]-=.005*6*np.pi*.001*a;return A,b,active
ref.system=system;q=ref.solve(x,x,.001);cpp=np.array(out['q']).reshape(2,6);ev=float(np.max(abs(q[:,:3]-cpp[:,:3])));ew=float(np.max(abs(q[:,3:]-cpp[:,3:])));assert out['wall_constraints']==2 and out['total_constraints']==3;assert ev<=1e-10 and ew<=1e-7
result={'status':'PASS','wall_constraints':2,'pair_constraints':1,'max_velocity_error_m_s':ev,'max_omega_error_rad_s':ew,'production_residual':out['residual'],'independent_method':'explicit KKT active-face enumeration, full coupled 12x12 resistance with bulk once'};(R/'validation/JOINT_CONSTRAINT_STATIC_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
