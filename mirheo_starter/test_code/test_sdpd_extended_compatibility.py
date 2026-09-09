"""Regression: historical CPU certificates must not prevent read-only reanalysis."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from py_scripts.fluid_physics.common import write_json
from py_scripts.sdpd_diagnostics.equilibration import cpu_validation


class HistoricalCertificate(unittest.TestCase):
    def test_stale_optional_certificate_is_not_current_evidence(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/'record.json'
            write_json(p, {'status':'PASS', 'code_sha256':{}})
            with patch('py_scripts.sdpd_diagnostics.equilibration.code_identity', return_value={'new.py':'digest'}):
                self.assertIsNone(cpu_validation({'cpu_validation_record':str(p)}))
                with self.assertRaisesRegex(ValueError, 'CPU_VALIDATION_CODE_SET_CHANGED'):
                    cpu_validation({'cpu_validation_record':str(p)}, required=True)
