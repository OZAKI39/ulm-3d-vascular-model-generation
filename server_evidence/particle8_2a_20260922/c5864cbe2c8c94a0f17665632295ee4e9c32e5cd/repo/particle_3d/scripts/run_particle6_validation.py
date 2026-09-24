#!/usr/bin/env python3
"""Generate actual LAMMPS CSV/JSON evidence and binary restart bundles."""
from pathlib import Path
import sys,argparse
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.particle6_validation import *
from particle_3d.particle3_cases import write_json,write_rows
from particle_3d.lammps_state import schema_contract
REPORT=PACKAGE/'reports/particle6';DATA=REPORT/'data'


def save(name,data):
    write_json(DATA/(name+'.json'),data)
    rows=data if isinstance(data,list) else [data]
    if isinstance(data,dict):
        if isinstance(data.get('errors'),list):rows=data['errors']
        elif 'exact_gap_records' in data:rows=data['exact_gap_records']
        elif 'rows' in data:rows=data['rows']
        elif 'before' in data and 'after' in data:
            rows=[dict(particle_id=a['particle_id'],before=a,after=b) for a,b in zip(data['before'],data['after'])]
        elif 'max_errors' in data:
            a=data['standalone'];b=data['bridge']
            rows=[dict(row=i,column=j,R_standalone=a['R'][i,j],R_bridge=b['R'][i,j],difference=a['R'][i,j]-b['R'][i,j]) for i in range(len(a['R'])) for j in range(len(a['R']))]
    write_rows(DATA/(name+'.csv'),rows)
    print(name,'saved',flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--part',choices=['all','synthetic','real','restart'],default='all');args=parser.parse_args()
    mu,viscosity=viscosity_from_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    provenance=frozen_provenance(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    if args.part in ['all','synthetic']:
        p,policy=mixed_particles()
        with LammpsParticleBridge(p[::-1],policy) as bridge:
            save('00_lammps_environment',environment(bridge))
            save('00_scope',dict(MODEL_FLAGS,viscosity=viscosity,provenance=provenance,schema=schema_contract()))
            save('01_state_roundtrip',dict(before=[x.to_dict() for x in p],after=[x.to_dict() for x in bridge.read()],errors=state_errors(p,bridge.read())))
            bridge.dump(REPORT/'logs/mixed_state.dump')
        save('02_neighbor_equivalence',neighbor_case(mu))
        sphere=run_parity(mu);mixed=run_parity(mu,True)
        save('03_one_step',dict(sphere=sphere['errors'][0],mixed=mixed['errors'][0],before=sphere['standalone'][0],
            standalone=sphere['standalone'][1],bridge=sphere['bridge'][1]))
        save('04_sphere_multistep',sphere);save('04_mixed_multistep',mixed)
        save('05_resistance',resistance_case(mu));save('10_neighbor_rebuild',rebuild_stress())
    if args.part in ['all','restart']:
        for mixed,name in [(True,'mixed'),(False,'sphere')]:
            result=restart_case(mu,REPORT/'checkpoints'/name,REPO,provenance,mixed=mixed)
            save('07_restart_'+name,result)
        # All three modes come from the mixed checkpoint.
        save('08_metadata_restart',json.loads((DATA/'07_restart_mixed.json').read_text())['metadata'])
    if args.part in ['all','real']:save('09_real_two_mb',real_case(REPO,mu))
    if all((DATA/(x+'.json')).exists() for x in ['04_sphere_multistep','04_mixed_multistep','07_restart_mixed','07_restart_sphere','09_real_two_mb','10_neighbor_rebuild']):
        rows=[];commands=[]
        for name in ['04_sphere_multistep','04_mixed_multistep','07_restart_mixed','07_restart_sphere','09_real_two_mb','10_neighbor_rebuild']:
            d=json.loads((DATA/(name+'.json')).read_text());rows.extend(dict(case=name,**row) for row in d['force_audits'])
            commands.extend(dict(case=name,command=c) for c in d['commands'])
        save('06_force_audit',rows);save('06_command_audit',commands)

if __name__=='__main__':main()
