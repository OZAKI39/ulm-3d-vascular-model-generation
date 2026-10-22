"""Fixed-window future runs and read-only historical population accounting."""
from pathlib import Path
from collections import Counter
import gzip, json, subprocess, sys
from .paths import ROOT, local_path, read_json, write_json, fresh_directory, sha256

def environment_spec(root=ROOT):
    spec=read_json(root/'configs/population/environment.json')
    for name in ('mesh','flow','sonovue_root'): spec[name]=str(local_path(spec[name],root))
    spec['boundaries']={k:str(local_path(v,root)) for k,v in spec['boundaries'].items()}
    return spec

def load_environment(root=ROOT):
    from particle_3d.population_v1.environment import load_environment as load
    return load(environment_spec(root))

def ledger(path):
    with gzip.open(path,'rt') as stream:
        for line in stream:
            if line.strip():yield json.loads(line)

def historical_summary(root=ROOT):
    run=root/'results/population/realization_0000_contactfix'
    summary=read_json(run/'SUMMARY.json')
    return dict(protocol='HISTORICAL_OUTCOME_DEPENDENT_COVERAGE_STOP',
                accepted=summary['accepted_birth_count'],outlets=summary['outlets'],
                censored=summary['censored_count'],numeric_failures=summary['solver_failed_count'],
                scientific_integrity_pass=summary['scientific_integrity_pass'],
                infusion_s=summary['infusion_stop_time_s'],end_s=summary['physical_time_s'])

def prepare_run(output, backend='cuda', root=ROOT):
    """Prepare only; output/results must not exist and remain distinct from frozen data."""
    from particle_3d.population_v1.fixed_window import validate_contract
    output=fresh_directory(output,root)
    c=read_json(root/'configs/population/FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json')
    validate_contract(c)
    write_json(output/'environment.json',environment_spec(root))
    write_json(output/'contract.json',c)
    source_identity={str(p.relative_to(root)):sha256(p) for p in sorted((root/'src').rglob('*')) if p.suffix in {'.py','.cpp','.cu'} and p.is_file()}
    write_json(output/'SOURCE_IDENTITY.json',source_identity)
    command=[sys.executable,'-m','particle_3d.population_v1.run_fixed_window_population',
             '--environment',str(output/'environment.json'),'--contract',str(output/'contract.json'),
             '--output',str(output/'results'),'--backend',backend,
             '--git-commit','independent-refactor-source-sha256:'+sha256(output/'SOURCE_IDENTITY.json')]
    write_json(output/'PLAN.json',dict(command=command,executed=False,protocol=c['schema'],historical_result_reused=False))
    return command

def run(command, output):
    with (Path(output)/'execution.log').open('x') as stream:
        return subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True).returncode
