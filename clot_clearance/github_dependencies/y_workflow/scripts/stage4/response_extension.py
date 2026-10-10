"""Offline coverage extension of the frozen APECSS model, only after measured shortfall."""
import json
import numpy as np
from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBuilder, AdaptiveResponseTableBubbleModel, ResponseNodeCache, response_errors
from scripts.stage4.common import ROOT, new_run, require, sha, write, finish


def main():
    require("acoustic")
    audit = json.loads((ROOT / "results/stage4/phase_screening/STATE.json").read_text())
    assert audit["status"] == "WARNING" and any(c["lumen_fraction_below_table"]>0 for c in audit["cases"].values())
    minimum = min(c["lumen_pressure_range_Pa"][0] for label,c in audit["cases"].items() if label != "B0")
    if minimum <= 0:
        raise RuntimeError("Exact acoustic pressure zero requires a separately validated zero-response endpoint; no artificial pressure floor")
    lower = float(10**(np.floor(np.log10(minimum))-1))
    bounds = [[1.779e-6,1.781e-6],[2.499e6,2.501e6],[lower,40000.]]
    config = dict(bounds=bounds,profile="production",workers=6,independent_holdouts=128,
                  reason="Measured Y acoustic domain contains pressures below the Stage2.5 5 kPa table minimum",
                  domain_audit_run=audit["run"],model="Unchanged APECSS/RPAR/Marmottant and parameter manifold; no fitted or extrapolated response")
    out = new_run("acoustic/response_extension",config,["acoustic"])
    builder = AdaptiveResponseTableBuilder(out,ROOT / "third_party/build/apecss/ultrasound_apecss",
        ROOT / "configs/stage2/sonovue_parameters.json",bounds=bounds,profile="production",workers=6,
        max_depth=24,max_nodes=12000,cache=ROOT / "results/stage2_5/response_table/node_cache")
    result = builder.build()
    model = AdaptiveResponseTableBubbleModel(out / "adaptive_response.h5")
    cache = ResponseNodeCache(builder.cache,builder.exe)
    rng = np.random.default_rng(20261006)
    points = rng.uniform(builder.bounds[:,0],builder.bounds[:,1],(128,3))
    points[:64,2] = np.exp(rng.uniform(np.log(lower),np.log(40000.),64))
    # Half the holdouts lie on the actual monodisperse/frequency operating line.
    points[::2,:2] = [1.78e-6,2.5e6]
    rows = []
    for point in points:
        node = cache.get(builder.case(point))
        try:
            assert node["status"] == "COMPLETE" and node["metrics"]["valid"] and node["metrics"]["harmonic_energy_resolved"]
            query = model.query(*point)
            errors = response_errors(query["harmonics_m3"],node,point[1],float(query["max_R_over_R0"]),float(query["min_R_over_R0"]))
            row = dict(point=point.tolist(),status="PASS" if max(errors.values())<=builder.config["tolerance"] else "FAIL",
                       errors={k:float(v) for k,v in errors.items()},cache_key=node["key"])
        except (AssertionError,ValueError) as error:
            row = dict(point=point.tolist(),status="FAIL",error=repr(error))
        rows.append(row)
        write(out / "HOLDOUTS.json",rows)
    passed = result["status"] == "PASS_EMPIRICAL_DOMAIN" and all(r["status"]=="PASS" for r in rows)
    if not passed:
        write(out / "TABLE_REJECTION.json",dict(reason="Incomplete adaptive coverage or independent holdout failure",build=result))
    finish(out,dict(status="PASS_WITH_SCOPE" if passed else "FAIL",builder=result,holdout_count=len(rows),
                   holdout_failures=sum(r["status"] != "PASS" for r in rows),bounds=bounds,
                   table=str((out / "adaptive_response.h5").relative_to(ROOT)),table_sha256=sha(out / "adaptive_response.h5"),
                   original_table_unchanged=True,scope="Narrow Stage4 operating manifold and new pressure coverage; empirical interpolation validation, no experimental recalibration"))


if __name__ == "__main__":
    main()
