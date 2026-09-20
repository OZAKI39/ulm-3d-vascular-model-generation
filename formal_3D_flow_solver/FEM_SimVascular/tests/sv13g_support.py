import json
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'reports/sv1_3g'
def load(name):return json.loads((R/(name+'.json')).read_text())
def actual(name):
    d=load(name)
    if d['status'] in ('NOT_RUN','BLOCKED'):pytest.skip(d.get('reason','Prerequisite not met'))
    assert d['status']=='PASS',d
    return d
