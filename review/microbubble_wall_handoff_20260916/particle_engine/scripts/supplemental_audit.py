from pathlib import Path
import json,re,hashlib,datetime,sys
R=Path(sys.argv[1]);ct=json.loads((R/'contracts/MINIMAL_LAMMPS_PARTICLE_ENGINE_CONTRACT.json').read_text());checks=[]
for sp in ct['cases']:
 p=R/'cases'/sp['case'];log=(p/'log.lammps').read_text();inp=(p/'in.lammps').read_text()
 # LAMMPS prints "Dangerous builds not checked" with check=no, not a zero counter.
 builds=[int(x) for x in re.findall(r'Neighbor list builds = (\d+)',log)]
 danger=[x.strip() for x in log.splitlines() if 'Dangerous builds' in x]
 rows=[x.split() for x in log.splitlines() if re.match(r'^\s*\d+\s+\d+\s+[\deE.+-]+\s+[\deE.+-]+\s*$',x)]
 thermo_steps=[int(x[0]) for x in rows];thermo_counts=[int(x[1]) for x in rows]
 check={'case':sp['case'],'requested_steps':sp['steps'],'neighbor_list_build_counts':builds,'actual_dangerous_build_lines':danger,'dangerous_build_count':'NOT_CHECKED_BY_LAMMPS','reason':'check no with every 1 forces unconditional rebuild; do not interpret absent integer count as measured zero','every_step_rebuild_configured':'neigh_modify every 1 delay 0 check no' in inp,'actual_rebuild_count_matches_steps':bool(builds) and builds[-1]==sp['steps'],'thermo_particle_counts':sorted(set(thermo_counts)),'thermo_count_constant':bool(thermo_counts) and all(n==sp['N'] for n in thermo_counts),'thermo_last_step_correct':bool(thermo_steps) and thermo_steps[-1]==sp['steps'],'status':'PASS'}
 if not all(check[x] for x in ['every_step_rebuild_configured','actual_rebuild_count_matches_steps','thermo_count_constant','thermo_last_step_correct']):check['status']='FAIL'
 checks.append(check)
result={'status':'PASS' if all(x['status']=='PASS' for x in checks) else 'FAIL','purpose':'Additional read-only audit of native neighbor diagnostics and all thermo particle counts; frozen finalizer and gates unchanged','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cases':checks}
(R/'validation/SUPPLEMENTAL_NATIVE_DIAGNOSTICS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'supplemental_audit':result['status']}));assert result['status']=='PASS'
