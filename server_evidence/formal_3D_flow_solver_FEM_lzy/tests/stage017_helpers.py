import json
from pathlib import Path
from copy import deepcopy
from functools import lru_cache
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPORT=ROOT/'reports/stage01_7';OUT=ROOT/'outputs/stage01_7'
def read(p):return json.loads(p.read_text())
POLICY=read(REPORT/'acceptance_policy.json');BASE=read(REPORT/'baseline_recomputed.json');CONTRACT=read(ROOT/'reports/stage01_6/planar_port_contract_v2.json')
def result():return read(OUT/'optimizer_result.json')
def source():return dict(np.load(ROOT/'inputs/stage01/tagged_surface_si.npz'))
def square():return np.array([[0.,0.,0.],[1,0,0],[1,1,0],[0,1,0],[.5,.5,0]]),np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]])
def measured_good():
 d=deepcopy(BASE);d['quality']['total_below_0_1']=0;d['quality']['cap_adjacent_below_0_1']=0
 d['quality']['low_quality_nearest_boundary_counts']={n:0 for n in BASE['quality']['low_quality_nearest_boundary_counts']}
 return d
@lru_cache(None)
def selected_recomputed():
 from fem3d.adaptive_qc import measure_volume
 return measure_volume(dict(np.load(OUT/'selected/mesh/volume_mesh.npz')),dict(np.load(OUT/'selected/surface/tagged_surface_si.npz')),source(),CONTRACT,POLICY)
