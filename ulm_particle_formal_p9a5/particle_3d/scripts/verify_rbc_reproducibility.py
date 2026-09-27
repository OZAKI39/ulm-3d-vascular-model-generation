#!/usr/bin/env python3
"""Two independent processes, identical explicit N/seed; keep receipt, not duplicates."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

PACKAGE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PACKAGE/"src"))
from particle_3d.rbc_distribution import digest


def verify(n,seed):
    records=[]
    with tempfile.TemporaryDirectory(prefix="particle2_repro_") as temp:
        for index in range(2):
            folder=Path(temp)/str(index)
            cmd=[sys.executable,"-B",str(PACKAGE/"scripts/generate_rbc_population.py"),
                 "--n",str(n),"--seed",str(seed),"--output",str(folder)]
            result=subprocess.run(cmd,check=True,text=True,capture_output=True)
            records.append(json.loads(result.stdout))
        csv=f"C57BL6_RBC_GEOMETRY_VALIDATION_{n}.csv"
        identical=(Path(temp)/"0"/csv).read_bytes()==(Path(temp)/"1"/csv).read_bytes()
        ledger_identical=(Path(temp)/"0/candidate_ledger.csv").read_bytes()==(Path(temp)/"1/candidate_ledger.csv").read_bytes()
    hashes=["diameter_array_sha256","volume_array_sha256","axes_array_sha256","sample_structured_array_sha256","csv_sha256"]
    passed=identical and ledger_identical and all(records[0][key]==records[1][key] for key in hashes)
    return dict(status="PASS" if passed else "FAIL",N=n,seed=seed,independent_process_count=2,
                csv_byte_identical=identical,candidate_ledger_identical=ledger_identical,
                float64_arrays_identical=all(records[0][key]==records[1][key] for key in hashes[:-1]),
                process_records=records,cross_numpy_version_byte_identity_claimed=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n",type=int,required=True);p.add_argument("--seed",type=int,required=True)
    p.add_argument("--output",type=Path,required=True);a=p.parse_args()
    record=verify(a.n,a.seed)
    a.output.write_text(json.dumps(record,ensure_ascii=False,indent=2)+"\n")
    print(record["status"],flush=True)
    if record["status"]!="PASS":raise SystemExit(1)


if __name__=="__main__":main()
