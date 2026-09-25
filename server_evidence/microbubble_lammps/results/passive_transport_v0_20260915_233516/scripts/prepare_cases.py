from pathlib import Path
import json,sys,os,math,hashlib
import numpy as np
P=Path(__file__).resolve().parents[1];ct=json.loads((P/'contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json').read_text());pop=json.loads((P/'provenance/POPULATION_IDENTITY.json').read_text());g2=json.loads((P/'contracts/CASE2_GEOMETRY_CONTRACT.json').read_text());g5=json.loads((P/'contracts/CASE5_GEOMETRY_CONTRACT.json').read_text());fc=json.loads((P/'contracts/FROZEN_FLOW_FIELD_CONTRACT.json').read_text())
cases=[]
def add(name,case,kind,cadv,idx=None,ranks=1,engine='cpu',phase='production'):
 key='case'+str(case);s=pop[key];ids=s['selected_source_ids'] if idx is None else [s['selected_source_ids'][idx]];ds=np.array(s['diameters_um']) if idx is None else np.array([s['diameters_um'][idx]])
 if case==5:
  ids=[q['source_bubble_id'] for q in g5['placements']];ds=np.array([q['diameter_m']*1e6 for q in g5['placements']])
 N=len(ids);dt=cadv*ct['field_speed_statistics'][kind]['dt_base_s'];coarse=ct['case_duration_policy'][{'uniform':'uniform_coarse_steps','linear':'linear_coarse_steps','real':'real_coarse_steps'}[kind]];factor=round(.5/cadv);steps=coarse*factor
 if engine=='gpu':steps=ct['case_duration_policy']['gpu_steps']
 stride=max(2,math.ceil(coarse/1000))*factor if engine!='gpu' else 2
 if kind=='real':
  xyz=np.array([g2['spawn_m']]) if case==2 else np.array([q['position_m'] for q in g5['placements']]);lo=np.array(fc['origin_m']);hi=lo+(np.array([fc['nx'],fc['ny'],fc['nz']])-1)*fc['dx_m']
 elif N>1:xyz=(np.array(np.unravel_index(np.arange(N),(10,10,10))).T+.5)*8e-6-40e-6;lo=np.full(3,-64e-6);hi=-lo
 else:xyz=np.zeros((N,3));lo=np.full(3,-64e-6);hi=-lo
 field='fields/'+{'real':'FROZEN_FLOW_FIELD_V0.h5','uniform':'UNIFORM_FIELD.h5','linear':'LINEAR_FIELD.h5'}[kind];sp={'name':name,'case':case,'phase':phase,'c_adv':cadv,'field_kind':kind,'field':field,'field_sha256':hashlib.sha256((P/field).read_bytes()).hexdigest(),'population':key,'source_bubble_ids':ids,'N':N,'diameters_um':ds.tolist(),'initial_positions_m':xyz.tolist(),'box_min_m':lo.tolist(),'box_max_m':hi.tolist(),'steps':steps,'dt_s':dt,'max_time_s':steps*dt,'trajectory_stride':stride,'mpi_ranks':ranks,'engine':engine,'real_geometry':kind=='real','check_pair_overlap':case==5,'margin_m':3*ct['dx_m'] if kind=='real' else 0}
 p=P/'cases'/name;p.mkdir();(p/'CASE_CONTRACT.json').write_text(json.dumps(sp,indent=2)+'\n')
 with (p/'particles.data').open('w') as f:
  f.write(f'Passive overdamped sphere storage only; density not used by integrator\n\n{N} atoms\n1 atom types\n\n')
  for a,k in enumerate('xyz'):f.write(f'{lo[a]:.17g} {hi[a]:.17g} {k}lo {k}hi\n')
  f.write('\nAtoms # sphere\n\n');np.savetxt(f,np.column_stack((np.arange(1,N+1),np.ones(N),ds*1e-6,np.ones(N),xyz)),fmt=['%d','%d']+['%.17g']*5)
  f.write('\nVelocities\n\n');np.savetxt(f,np.column_stack((np.arange(1,N+1),np.zeros((N,6)))),fmt=['%d']+['%.17g']*6)
 rel=os.path.relpath(P/field,p);wall='../../provenance/closed_geometry_m.stl' if kind=='real' else 'NONE'
 lines=['units si','atom_style sphere','boundary p p p','newton off',('processors 2 2 1' if ranks==4 else 'processors 1 1 1'),'read_data particles.data','comm_modify vel yes','pair_style zero 6e-6','pair_coeff * *','neighbor 1e-6 bin','neigh_modify every 1 delay 0 check no',f'fix passive all sonovue/passive/flow {rel} {stride} trajectory {wall} {sp["margin_m"]:.17g} {int(case==5)}',f'timestep {dt:.17g}',f'thermo {max(1,steps//10)}','thermo_style custom step atoms time','thermo_modify lost error norm no format float %.17g','run 0','write_dump all custom initial.dump id diameter xu yu zu vx vy vz modify sort id format float %.17g',f'run {steps}','write_dump all custom final.dump id diameter xu yu zu vx vy vz modify sort id format float %.17g','print "PASSIVE_CASE_COMPLETED"']
 (p/'in.lammps').write_text('\n'.join(lines)+'\n');cases.append(sp)
for label,cadv in [('c050',.5),('c025',.25),('c0125',.125)]:
 add('case0_'+label,0,'uniform',cadv,phase='timestep_validation');add('case1_'+label,1,'linear',cadv,phase='timestep_validation')
 if g2['status']=='PASS':add('case2_'+label,2,'real',cadv,phase='timestep_validation')
for label,cadv in [('c025',.25),('c0125',.125)]:
 for i,q in enumerate(['d10','d50','d90']):add('case3_'+q+'_'+label,3,'linear',cadv,idx=i)
 for n in [1,4]:add(f'case4_mpi{n}_'+label,4,'uniform',cadv,ranks=n)
 if g5['N']>=3:add('case5_'+label,5,'real',cadv)
 add('kokkos_'+label,4,'uniform',cadv,engine='gpu')
(P/'contracts/CASE_INDEX.json').write_text(json.dumps(cases,indent=2)+'\n')
print('CASE_INPUTS_PREPARED',len(cases))
