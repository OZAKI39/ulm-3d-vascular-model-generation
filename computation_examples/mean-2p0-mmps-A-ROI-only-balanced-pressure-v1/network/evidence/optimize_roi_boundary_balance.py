"""ROI-only design in strict regression -> basis -> cross-check -> freeze order."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from network_1d0d.audit import write_json,write_csv
from network_1d0d.roi_only_hydraulics import load_roi_cache
from network_1d0d.roi_boundary_design import regress_h0,design,freeze_design


def run(output,freeze=False):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    h0=ROOT/'reports/a_network_1d0d_boundary_v2_idealized/data'
    cache=load_roi_cache(h0/'analysis_A_H0_graph_si.npz',ROOT/'reports/a_network_1d0d_boundary_v1/data/roi_ports_in_a.json')
    write_json(output/'roi_geometry_cache_summary.json',cache.summary())
    np.savez_compressed(output/'roi_geometry_cache.npz',**{k:v for k,v in vars(cache).items() if isinstance(v,np.ndarray)})
    regression=regress_h0(cache,json.loads((h0/'roi_operating_point_H0.json').read_text()))
    write_json(output/'roi_mixed_bc_regression.json',regression)
    cert=json.loads((ROOT/'reports/balanced_three_outlet_forward_v1/feasibility_certificate.json').read_text())
    solution,basis,trace=design(cache,regression,cert)
    write_json(output/'roi_pressure_basis.json',basis);write_json(output/'roi_equal_split_solution.json',solution)
    write_csv(output/'roi_iteration_trace.csv',trace)
    if freeze:freeze_design(output/'roi_boundary_design_frozen.yaml',solution)
    return solution


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--freeze',action='store_true')
    a=p.parse_args();print(json.dumps(run(a.output,a.freeze),indent=2))
