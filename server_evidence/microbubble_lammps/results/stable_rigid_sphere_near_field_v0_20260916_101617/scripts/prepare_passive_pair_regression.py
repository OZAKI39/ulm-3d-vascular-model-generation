from pathlib import Path
import csv,json,hashlib
R=Path(__file__).resolve().parents[1];old=Path('/home/lzy/projects/compre_output/passive_transport_v0/20260915_233516/cases/case4_mpi1_c025');orig=json.loads((old/'CASE_CONTRACT.json').read_text());D=R/'cases/J2_passive_separated_pair';D.mkdir(exist_ok=True)
c=json.loads((R/'cases/J_passive/CASE_CONTRACT.json').read_text());c.update(name=D.name,N=2,radii_m=[x*.5e-6 for x in orig['diameters_um'][:2]],initial_positions_m=orig['initial_positions_m'][:2],ids=[1,2],field=orig['field'],field_sha256=orig['field_sha256'],dt_max=orig['dt_s'],max_time=orig['max_time_s'],wall=False,source='original passive case4 MPI1, selected original IDs 1/2; both remain inside original box, no wrapping; nonempty sphere pair outside lubrication cutoff')
(D/'CASE_CONTRACT.json').write_text(json.dumps(c,indent=2)+'\n');data='Exact original passive case4 IDs 1/2; overdamped storage\n\n2 atoms\n1 atom types\n\n'
for axis in 'xyz':data+=f'-6.4e-05 6.4e-05 {axis}lo {axis}hi\n'
data+='\nAtoms # sphere\n\n'
for tag,a,x in zip(c['ids'],c['radii_m'],c['initial_positions_m']):data+=f'{tag} 1 {2*a:.17g} 1 '+' '.join(format(v,'.17g') for v in x)+'\n'
data+='\nVelocities\n\n1 0 0 0 0 0 0\n2 0 0 0 0 0 0\n';(D/'particles.data').write_text(data)
values={'field':'../../'+c['field'],'data':'particles.data','dt_max':c['dt_max'],'c_gap':.4,'max_time':c['max_time'],'margin':0,'stride':1,'drive_mode':0,'drive_scale':0,'steps_cap':200000,'wall':'NONE'};(D/'case.cfg').write_text(''.join(f'{k} {v}\n' for k,v in values.items()))
idx=json.loads((R/'contracts/CASE_INDEX.json').read_text());assert not any(x['name']==D.name for x in idx);idx.append(c);(R/'contracts/CASE_INDEX.json').write_text(json.dumps(idx,indent=2)+'\n')
source=old/'trajectory_rank0.csv';out=R/'provenance/J2_OLD_SELECTED_TRAJECTORIES.csv';count=0
with open(source) as f,open(out,'w') as g:
 reader=csv.reader(f);writer=csv.writer(g);header=next(reader);writer.writerow(header)
 for row in reader:
  if int(row[0]) in [1,2]:writer.writerow(row);count+=1
meta={'role':'NONEMPTY_PAIR_OUTSIDE_CUTOFF_REGRESSION','source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'original_ids':[1,2],'rows':count,'selected_csv_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'reason':'J single sphere is retained, but cannot alone demonstrate a nonempty pair outside the cutoff; supplement completes requested coverage','physics_contract_changed':False,'source_changed':False}
(R/'contracts/J2_REGRESSION_COVERAGE.json').write_text(json.dumps(meta,indent=2)+'\n');print('J2_PREPARED',count)
