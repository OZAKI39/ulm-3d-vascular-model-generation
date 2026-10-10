"""Final regression and the original immutable Stage 3 evidence manifest."""
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import h5py
import numpy as np
import yaml
from scripts.stage4.common import ROOT,new_run,require,sha,write,finish


def stage3_metadata(commit,out):
    names = ["PHYSICS_ACCEPTANCE.json","docs/STAGE3_VALIDATION_REPORT.md",
             "docs/HISTORY_FORCE_SENSITIVITY.md","docs/KNOWN_LIMITATIONS.md",
             "docs/METHODS_MODEL_DESCRIPTION.md"]
    rows = []
    for name in names:
        original = subprocess.check_output(["git","show",f"{commit}:{name}"],cwd=ROOT)
        expected = hashlib.sha256(original).hexdigest()
        actual = sha(ROOT/name) if (ROOT/name).is_file() else None
        rows.append(dict(path=name,expected_sha256=expected,actual_sha256=actual,unchanged=actual==expected))
    result = dict(status="PASS" if all(row["unchanged"] for row in rows) else "FAIL",
                  baseline_commit=commit,files=rows)
    write(out/"STAGE3_BASELINE_METADATA.json",result)
    return result


def stage4_identities(out):
    flow,mesh,acoustic,controls = [require(name) for name in ["simvascular/steady","simvascular/mesh","acoustic","phase_screening"]]
    paths = [row["run"] for name in ["primary_only","statistics/paired_runs","primary_secondary"] for row in require(name)["rows"]]
    paths += [str(path.parent.relative_to(ROOT)) for path in (ROOT/require("convergence")["run"]).rglob("RUN_IDENTITY.json")]
    expected = {"flow_sha256":flow["snapshot_sha256"],"mesh_sha256":mesh["mesh_sha256"],
                "acoustic_basis_sha256":acoustic["basis_sha256"],"response_table_sha256":sha(ROOT/controls["response_table"]),
                "geometry_hash":flow["geometry_hash"]}
    files = {flow["snapshot"]:expected["flow_sha256"],mesh["mesh"]:expected["mesh_sha256"],
             acoustic["basis"]:expected["acoustic_basis_sha256"],controls["response_table"]:expected["response_table_sha256"]}
    production = {str(p.relative_to(ROOT)):sha(p) for p in (ROOT/"src").rglob("*")
                  if p.is_file() and p.suffix in {".py",".yaml"} and "__pycache__" not in p.parts}
    gold=ROOT/"results/stage4/performance_refactor/gold_baseline"
    gold_source=json.loads((gold/"cases/B2/SOURCE_IDENTITY.json").read_text())["source_sha256"]
    original_physics={name:digest for name,digest in gold_source.items() if name.startswith("src/")}
    assert all(production.get(name)==digest for name,digest in original_physics.items()),"Original production physics changed"
    decision=require("performance_refactor/backend_decision")
    manifest=json.loads((gold/"GOLD_BASELINE_MANIFEST.json").read_text())
    assert sha(gold/"GOLD_BASELINE_MANIFEST.json")==decision["gold_manifest_sha256"]
    gold_changed=[name for name,row in manifest["files"].items()
        if sha(gold/name)!=row["sha256"] or (gold/name).stat().st_size!=row["bytes"]]
    gold_writable=[str(path.relative_to(gold)) for path in gold.rglob("*") if path.stat().st_mode&0o222]
    assert not gold_changed and not gold_writable,(gold_changed,gold_writable)
    write(out/"GOLD_BASELINE_PRESERVATION.json",dict(status="PASS",verified_files=len(manifest["files"]),
        manifest_sha256=sha(gold/"GOLD_BASELINE_MANIFEST.json"),changed=gold_changed,writable_paths=gold_writable))
    checked = []
    for path in sorted(set(paths)):
        folder = ROOT/path
        assert not (folder/"TERMINATION_POLICY_VIOLATION.json").exists(),(path,"Objective termination violation")
        identity = json.loads((folder/"RUN_IDENTITY.json").read_text())
        assert all(identity[key]==value for key,value in expected.items()),path
        source = json.loads((folder/"SOURCE_IDENTITY.json").read_text())["source_sha256"]
        recorded_production = {name:digest for name,digest in source.items() if name.startswith("src/")}
        assert all(recorded_production.get(name)==digest for name,digest in original_physics.items()),(path,"Original physics source differs from gold")
        assert all(production.get(name)==digest for name,digest in recorded_production.items()),(path,"Versioned backend source changed after its formal run")
        certificate = identity["config"].get("performance_execution")
        if certificate:
            if "wall_projection" in certificate:
                helper = "scripts/stage4/geometry_certificate.py"
                assert certificate["wall_projection"] == "per_particle_interior_certificate_v1"
                assert certificate["helper_sha256"] == source[helper] == sha(ROOT/helper)
                proof = json.loads((ROOT/certificate["equivalence_evidence"]/"SUMMARY.json").read_text())
                assert proof["status"] == "PASS"
            if "backend" in certificate:
                assert certificate["backend"]==decision["recommended_backend_for_JK"]
                assert certificate["backend_decision_sha256"]==sha(ROOT/decision["run"]/"SUMMARY.json")
        with np.load(folder/"INITIAL_POPULATION.npz") as initial:
            actual = hashlib.sha256(initial["positions_m"].tobytes()+initial["radii_m"].tobytes()+initial["particle_id"].tobytes()).hexdigest()
            assert actual==identity["initial_population_sha256"],path
            count = len(initial["particle_id"])
        branch = json.loads((folder/"BRANCH_METRICS.json").read_text())
        assert sum(branch[k] for k in ["N_A","N_B","N_unassigned","N_wall_loss"])==count
        with h5py.File(folder/"particles/trajectories.h5") as h:
            assert h.attrs["complete"] and len(h["particle_id"])==count
            assert np.isfinite(h["positions_m"][:]).all()
            assert float(h["time_s"][-1])>=identity["config"]["particles"]["duration_s"]-1e-12 or not h["active"][-1].any()
        if "termination_policy" in identity["config"]:
            from mevo_mb.coupling.termination import SimulationTerminationPolicy
            completed=json.loads((folder/"particles/SUMMARY.json").read_text())
            policy=SimulationTerminationPolicy(**identity["config"]["termination_policy"])
            reason=policy.reason(time_s=completed["duration_s"],
                observation_limit_s=identity["config"]["particles"]["duration_s"],
                total_particles=count,terminally_assigned=branch["N_assigned"])
            assert reason is not None and completed["termination_reason"]==reason,(path,"Termination policy not satisfied")
        if not identity["acoustic_OFF"]:
            field_path = identity.get("acoustic_field_path") or controls["cases"][identity["control_case"]]["field"]
            files[field_path] = identity["acoustic_field_sha256"]
        checked.append(dict(run=path,initial_hash_verified=actual,bookkeeping_complete=True,
                            production_source_sha256_verified=True,
                            wall_certificate_helper_sha256=certificate.get("helper_sha256") if certificate else None,
                            original_physics_source_files=len(original_physics),recorded_source_files=len(recorded_production),
                            versioned_performance_backend=certificate.get("backend","original") if certificate else "original"))
    mismatched = [path for path,digest in files.items() if sha(ROOT/path)!=digest]
    assert not mismatched,mismatched
    config = yaml.safe_load((ROOT/"configs/stage4/bifurcation.yaml").read_text())
    workflow = yaml.safe_load((ROOT/"configs/stage4/workflow.yaml").read_text())
    sys.path.insert(0,str(Path(workflow["solver_executable"]).parents[1]/"vendor"))
    from flow_solver_support.transient_checks import checkpoint_audit
    final = ROOT/flow["snapshot"]
    step = int(final.stem.rsplit("_",1)[1])
    native = checkpoint_audit(final.with_name(f"stFile_{step:03d}.bin"),step,
                               config["flow"]["steady_relaxation_dt_s"],ranks=config["flow"]["mpi_ranks"])
    result = dict(status="PASS",particle_runs=checked,verified_input_sha256=files,
                  verified_production_source_sha256=production,unchanged_gold_physics_source_sha256=original_physics,
                  final_native_checkpoint=native,backend_decision=decision["run"])
    write(out/"STAGE4_RESULT_IDENTITIES.json",result)
    return result


def main():
    baseline = require("audit")
    out = new_run("audit/final",dict(baseline_run=baseline["run"]))
    protected = json.loads((ROOT / baseline["run"] / "STAGE3_BASELINE_FILES.json").read_text())
    changed = [name for name,digest in protected.items() if not (ROOT/name).is_file() or sha(ROOT/name)!=digest]
    assert subprocess.check_output(["git","rev-parse","stage3-minimal-validated^{}"],text=True).strip()==baseline["baseline_commit"]
    metadata = stage3_metadata(baseline["baseline_commit"],out)
    changed += [row["path"] for row in metadata["files"] if not row["unchanged"]]
    checks = []
    for name,command in [("regression",[".venv/bin/pytest","-q","--junitxml="+str(out/"pytest.xml")]),
                         ("archive",[".venv/bin/pytest","-q","-m","deprecated_reference"]),
                         ("ruff",[".venv/bin/ruff","check","src","tests","scripts","validation/history_force_sensitivity","validation/history_branch_robustness"])]:
        result = subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
        (out / (name+".log")).write_text(result.stdout+result.stderr)
        checks.append(dict(name=name,returncode=result.returncode,output=(result.stdout+result.stderr)[-1200:]))
    forbidden = []
    for path in (ROOT/"src").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            names = [node.module or ""] if isinstance(node,ast.ImportFrom) else [x.name for x in node.names] if isinstance(node,ast.Import) else []
            if any("history_branch_robustness" in name or "history_force_sensitivity" in name for name in names):
                forbidden.append(str(path.relative_to(ROOT)))
    identities = stage4_identities(out)
    finish(out,dict(status="PASS" if not changed and not forbidden and all(c["returncode"]==0 for c in checks) else "FAIL",
                    protected_files=len(protected),changed=changed,production_history_imports=forbidden,checks=checks,
                    stage3_acceptance_sha256=sha(ROOT/"PHYSICS_ACCEPTANCE.json"),
                    baseline_metadata_files=len(metadata["files"]),baseline_metadata_identity_check=metadata["status"],
                    stage4_completed_run_count=len(identities["particle_runs"]),stage4_input_identity_check="PASS",
                    stage4_production_source_count=len(identities["verified_production_source_sha256"]),
                    stage4_unchanged_gold_physics_source_count=len(identities["unchanged_gold_physics_source_sha256"]),
                    baseline_commit=baseline["baseline_commit"],scope="Regression and evidence preservation; numerical tests do not establish physical calibration"))


if __name__=="__main__":
    main()
