import json
import pytest
from sv_validation.sv12 import CONFIG,load,validate_frozen
import sv_validation.sv12 as sv12
from sv_validation.validation import ValidationError

def test_exact_options_inherited():
    config=json.loads((CONFIG/'petsc_options.json').read_text());old=load('petsc_short_execution','sv1_1')
    assert config['PETSC_OPTIONS']==old['PETSC_OPTIONS']
    assert config['mpi_ranks']==old['mpi_ranks']==4
    assert config['OMP_NUM_THREADS']==old['omp_num_threads']==1

def test_changed_petsc_options_rejected(tmp_path,monkeypatch):
    config=json.loads((CONFIG/'petsc_options.json').read_text())
    config['PETSC_OPTIONS']=config['PETSC_OPTIONS'].replace('-sub_pc_factor_levels 2','-sub_pc_factor_levels 3')
    (tmp_path/'petsc_options.json').write_text(json.dumps(config))
    monkeypatch.setattr(sv12,'CONFIG',tmp_path)
    with pytest.raises(ValidationError,match='FROZEN_INPUT_CHANGED: PETSc'):validate_frozen()

def test_early_mass_solver_ranking_removed_only_in_current_policy():
    assert json.loads((CONFIG/'policy.json').read_text())['early_transient_mass_solver_ranking'] is False
    assert load('stage_result','sv1_1')['reason']=='SHORT_RUN_MASS_NOT_IMPROVED'
