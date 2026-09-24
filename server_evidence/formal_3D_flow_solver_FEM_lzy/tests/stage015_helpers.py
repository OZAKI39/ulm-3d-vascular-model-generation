import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/'reports/stage01_5'
OUT=ROOT/'outputs/stage01_5'
CANDIDATES=['candidate_A','candidate_B','candidate_C']
def read(path): return json.loads(path.read_text())
def square():
    return np.array([[0.,0.,0.],[1,0,0],[1,1,0],[0,1,0],[.5,.5,0]]),np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]])
