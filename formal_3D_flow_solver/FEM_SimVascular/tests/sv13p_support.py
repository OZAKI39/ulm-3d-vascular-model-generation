import hashlib,json,math,re
from pathlib import Path
import pytest
from sv_validation.sv13p import *
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'reports/sv1_3p'
def read(name):return json.loads((R/(name+'.json')).read_text())
def accepted(name):
 d=read(name);assert d['status']=='PASS';return d
def cases():return [json.loads(p.read_text()) for p in R.glob('*_acceptance.json') if 'before_parser' not in p.name]
def candidate(name):return read('candidate_summary')['candidates'][name]
