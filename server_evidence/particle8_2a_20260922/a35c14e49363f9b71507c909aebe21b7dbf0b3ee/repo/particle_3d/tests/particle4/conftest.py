from pathlib import Path
import sys
import pytest
PACKAGE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PACKAGE/'src'))

@pytest.fixture(scope='session')
def p4_repo():return PACKAGE.parent
