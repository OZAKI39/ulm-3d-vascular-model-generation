from pathlib import Path
import pytest
from sv_validation.sv12 import ROOT,OUTPUT
from sv_validation.validation import parse_result_vtu,require,ValidationError
from sv12_support import artifact

def test_all_scheduled_states_present():
    execution=artifact('flow_execution');states=artifact('saved_state_qc')['states']
    assert [s['step'] for s in states]==list(range(10,execution['last_step']+1,10))
    assert execution['last_step'] in (400,800)
    for state in states:assert (ROOT/state['path']).is_file()

def test_missing_vtu_rejected(tmp_path):
    with pytest.raises((FileNotFoundError,OSError,ValueError)):parse_result_vtu(tmp_path/'missing.vtu')
