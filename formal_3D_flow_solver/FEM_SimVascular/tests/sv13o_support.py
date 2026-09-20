import hashlib,json,re,xml.etree.ElementTree as ET
from pathlib import Path
import pytest
from sv_validation.sv13o import *
from sv_validation.sv13 import stop_gate
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'reports/sv1_3o';C=ROOT/'configs/sv1_3o'
def read(name):return json.loads((R/(name+'.json')).read_text())
def accepted(name):
 d=read(name);assert d['status']=='PASS',d.get('errors',d);return d
def execution(name):return read('remote/'+name+'_execution')
def result(name):return read(name+'_acceptance')

