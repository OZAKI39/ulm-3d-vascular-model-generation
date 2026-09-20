from sv_validation.sv13 import check_reference
from sv_validation.provenance import sha256
from sv_validation.validation import ValidationError
from sv13_support import artifact
import pytest
def test_reference_content_is_unchanged():
    assert check_reference()
    assert artifact('reference_freeze')['designation']=='VALIDATION_REFERENCE'
def test_modified_reference_rejected(tmp_path):
    p=tmp_path/'vtu';p.write_text('original');expected=sha256(p);p.write_text('changed')
    with pytest.raises(ValidationError,match='REFERENCE_CHANGED'):
        check_reference([{'path':str(p),'sha256':expected}])
