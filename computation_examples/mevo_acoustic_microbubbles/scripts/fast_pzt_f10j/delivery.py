"""Generate layered decisions and reports directly from the measured evidence."""

import subprocess
import numpy as np
from .common import ROOT, BASE, BASELINE, js, sha, write


def main():
    audit = js(BASE / "analysis/ORIGINAL_AUDIT.json")
    valid = js(BASE / "optimized_v2/RESULTS.json")
    corrected = js(BASE / "corrected_v1/RESULTS.json")
    guards = js(BASE / "guards_v1/RESULTS.json")
    perf = js(BASE / "performance_v1/RESULTS.json")
    cold = js(BASE / "process_v1/RESULTS.json")
    row = valid["rows"][-1]
    new = row["gpu"]
    m = row["metrics"]
    original = ROOT / "results/fast_pzt/f10i"
    manifest = js(original / "ARTIFACT_MANIFEST.json")["files"]
    for name, entry in manifest.items():
        assert sha(original / name) == entry["sha256"], name
    inputs = js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json")
    for name, entry in inputs.items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == entry["sha256"], name
    # Use the explicit original worktree; the current worktree is independently based on this HEAD.
    history_root = type(ROOT)(
        "/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10i_20261012/computation_examples/mevo_acoustic_microbubbles"
    )
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=history_root, text=True).strip()
    status = subprocess.check_output(["git", "status", "--porcelain"], cwd=history_root, text=True)
    assert head == BASELINE and not status
    write(
        BASE / "BASELINE_AUDIT.json",
        dict(
            status="PASS",
            F10I_complete_HEAD=head,
            F10I_status_porcelain=status,
            F10I_artifacts_verified=len(manifest),
            F10I_manifest_sha256=sha(original / "ARTIFACT_MANIFEST.json"),
            old_v1_execution=js(original / "gpu_v1/IDENTITY.json")["execution_commit"],
            old_v2_execution=js(original / "gpu_v2/IDENTITY.json")["execution_commit"],
            F10I_failure_sha256=sha(original / "analysis/GPU_FAILURE_ANALYSIS.json"),
            F10I_gate_sha256=sha(original / "analysis/FINAL_GATE.json"),
            input_files_verified=len(inputs),
            worktrees=subprocess.check_output(["git", "worktree", "list"], cwd=ROOT, text=True),
            local_AGENTS="No applicable AGENTS.md in worktree or ancestors; prior projects not restarted",
            remote_AGENTS="/workspace/AGENTS.md -> /etc/vast-agents-guide.md; read before remote actions",
            server_guide_sha256=js(BASE / "ENVIRONMENT.json")["guide_sha256"],
        ),
    )
    names = ["original", "corrected", "optimized", "cpu"]
    medians = {}
    for n in names:
        samples = [next(v for v in s["measurements"] if v["name"] == n) for s in perf["samples"]]
        medians[n] = dict(
            median_wall_ms=float(np.median([v["matched_full_wall_s"] for v in samples]) * 1e3),
            min_wall_ms=min(v["matched_full_wall_s"] for v in samples) * 1e3,
            max_wall_ms=max(v["matched_full_wall_s"] for v in samples) * 1e3,
            median_preparation_solve_ms=float(
                np.median([v["solve_and_preparation_s"] for v in samples]) * 1e3
            ),
            median_quality_ms=float(np.median([v["quality_wall_s"] for v in samples]) * 1e3),
            median_output_ms=float(np.median([v["matched_output_commit_wall_s"] for v in samples]) * 1e3),
        )
    process_cost = {
        p.stem.replace("_PROCESS", ""): js(p)["full_process_wall_s"] for p in BASE.glob("*_PROCESS.json")
    }
    performance = dict(
        medians=medians,
        new_gpu_over_cpu_cost=medians["optimized"]["median_wall_ms"] / medians["cpu"]["median_wall_ms"],
        performance_pass=medians["optimized"]["median_wall_ms"] < medians["cpu"]["median_wall_ms"],
        scalar_reads={"original": 35, "corrected": 17, "optimized": 12},
        final_array_reads={"original": 12, "corrected": 14, "optimized": 14},
        final_read_count_method="Count of explicit numpy() calls in inherited solve plus two final F10J linear-audit calls; intermediate count is runtime measured.",
        independent_processes=cold,
        science_process_wall_s=process_cost,
        total_registered_process_wall_s=sum(process_cost.values()),
        pure_kernel_time_claim=False,
        production_speedup_claim=False,
    )
    write(BASE / "analysis/PERFORMANCE_SUMMARY.json", performance)
    assert not performance["performance_pass"]
    assert corrected["numerical_admitted"] and valid["numerical_admitted"] and guards["status"] == "PASS"
    gate = dict(
        mathematical_convergence="GPU_FIXED_STATE_NUMERICAL_PASS",
        GPU_implementation="FIXED_STATE_REFERENCE_AND_GUARD_PASS",
        small_fixture_tests_passed=sum(r["test_pass"] for r in valid["rows"] if r["n"] < 32),
        small_accepted_states=16,
        small_expected_quality_rejections=2,
        real_32_state_pass=True,
        nonlinear_gradient_bound_m=new["conservative_equation_bound_m"],
        original_gradient_target_m=1e-14,
        old_and_new_PCG_targets_unchanged=True,
        all_linear_true_and_recursive_checks_pass=True,
        position_and_force_metrics=m,
        invalid_input_and_budget_checks="PASS",
        incomplete_linear_solve_guard="PASS",
        performance="NO_B1_COST_ADVANTAGE",
        GPU_B2_PRODUCTION_NOT_ADMITTED=True,
        decision="STOP_GPU_B2_PRODUCTION_INTEGRATION; retain opt-in corrected module and CPU numerical reference",
        contact_physics="PHYSICS_UNCALIBRATED",
        stiffness_N_per_m=0.04641319518378658,
        group_morphology="GROUP_LEVEL_NUMERICAL_EVIDENCE_RETAINED_FROM_F10H_I",
        strict_instantaneous_edges="RAW_EDGE_TIMING_FAIL_RETAINED_12_OF_29_MAX_130_US",
        persistent_topology="PERSISTENT TOPOLOGY VALIDATION NOT RUN",
        historical_candidate="HISTORICAL_INITIAL_STATE_CANDIDATE_ONLY_F5_417P5_MS_NOT_B2_TRUTH",
        training_ready=False,
        production=False,
        labels=["EXPLORATORY", "PHYSICS_UNCALIBRATED"],
        new_trajectories=0,
        new_CFD=0,
        new_acoustic_fields=0,
        new_AI_training=0,
        execution_commits={p.parent.name: js(p)["execution_commit"] for p in BASE.glob("*/IDENTITY.json")},
        remote_archive=js(BASE / "REMOTE_VERIFICATION.json"),
    )
    write(BASE / "analysis/FINAL_GATE.json", gate)
    write(
        BASE / "analysis/ATTEMPT_LEDGER.json",
        dict(
            scientific_attempts=[
                dict(
                    label=k,
                    exit_code=js(BASE / (k + "_PROCESS.json"))["exit_code"],
                    kind="expected failure reproduction"
                    if k == "audit_v1"
                    else "negative guards included"
                    if k == "guards_v1"
                    else "completed",
                )
                for k in process_cost
            ],
            local_engineering_corrections=[
                dict(
                    step="read-only Warp inspection",
                    error="ImportError: _LinearSolverState; actual public name LinearSolverState",
                    effect="Inspection only, corrected before solver instrumentation",
                ),
                dict(
                    step="first offline audit write",
                    error="FileNotFoundError: missing analysis directory",
                    effect="Created directory and reran unchanged inputs; no scientific execution failed",
                ),
                dict(
                    step="offline trace labelling review",
                    error="Before-direction outer evaluation initially selected instead of immediately pre-trial evaluation",
                    effect="Corrected analysis before candidate preregistration; raw read-only GPU audit unchanged",
                ),
                dict(
                    step="historical schema v1",
                    error="KeyError: active",
                    log="topology_screen_v1_failure.log",
                    effect="Preserved log; explicit schema guard added",
                ),
                dict(
                    step="historical schema v2",
                    error="Incomplete coverage of newer trajectory/active schema",
                    log="topology_screen_v2.log",
                    effect="Mapped actual new schema; final v3 includes F10C-I",
                ),
            ],
            old_F10I_failures_preserved=True,
            no_failed_scientific_attempt_deleted=True,
        ),
    )
    docs = ROOT / "docs"
    linear_lines = "\n".join(
        f"|{r['kind']} / {r['outer']}|{r['recursive_final']:.6e}|{r['true_final']:.6e}|{r['tolerance']:.6e}|{'PASS' if r['true_pass'] and r['recursive_pass'] else 'FAIL'}|"
        for r in audit["linear"]
    )
    last = next(e for e in audit["evaluations"] if e["trial"] and e["outer"] == 4)
    audit_text = f"""# F10J GPU solver audit

The original saved 32-bubble failure is reproduced exactly. Both the base PCG and **the first contact-direction PCG** exceed their own targets. The last direction is accurate and descending; total-energy roundoff rejects all 20 Armijo candidates. The original mathematical stopping rule was retained and met after the correction.

## Actual equation and units

Let R=1.78e-6 m, dt=1e-5 s, zeta=6 pi mu R, mu=0.001 Pa s, and kc=0.04641319518378658 N/m (F10G 2kc, exploratory). The initial-step resistance is dimensionless A=I+sum c_ij B_old^T B_old. The code uses c=0.5[phi(max(h/R,0.001))-phi(0.1)] for h/R<0.1; phi(e)=0.5/e-0.45 log(e)-(3/56)e log(e). It is fixed during the contact solve. Free velocity already includes the archived external acoustic/background inputs; no new force or flow is generated.

With x_s=(x-origin)/R, b=dt*free/R, A u=b, pred=x_s+u, q dimensionless, alpha=dt kc/zeta, y=pred+q:

`E(q)=0.5 q^T A q + 0.5 alpha sum_pairs max(2-|y_j-y_i|,0)^2`.

`g=Aq-alpha f(y)`; physical contact force is `kc R f(y)`. The physical equation is `x_new-x_old = dt A^-1 [free + F_contact/zeta]`. The global-coordinate residual in metres is evaluated independently from actual returned positions and real pair geometry.

Every outer iteration updates pair distances, normals and the unilateral active set. A stays frozen. The direction solves `H d=-g`, `H=A+alpha B_active^T B_active`, using a block Jacobi inverse of each diagonal 3x3 block. H is a positive Gauss-Newton metric, **not** the full nonlinear Hessian and not a proof of global convexity or unique physical equilibrium. Measured symmetry errors are at most {max(r["symmetry_max"] for r in audit["linear"]):.3e}; min eigenvalues are approximately 1; all block inverses are positive definite; no used direction has nonpositive PCG curvature. A's maximum eigenvalue is {audit["linear"][0]["eigen_max"]:.6f}. CPU and GPU A differ by {audit["gpu_cpu_operator_max_diff"]:.3e} in absolute matrix entries from coordinate evaluation roundoff; this is not a physical model change.

## Three different residuals

1. PCG recursive r is updated by the recurrence; raw norms are dimensionless and compared with `max(1e-14,1e-12 ||b||_2)` for that solve's RHS.
2. True linear r is recomputed on the GPU as b-Ax or -g-Hd using the actual operator and GPU solution. The preconditioned diagnostic is sqrt(r^T M r), **not** the unpreconditioned stopping norm. Relative residual divides by RHS norm. All arrays and independent dense-reference linear errors are archived.
3. Full nonlinear position residual is independently recomputed in metres from the actual final contact forces and A inverse. It is not a PCG residual.

For exact base prediction, `R ||A^-1 g||_2 <= R ||g||_2` since lambda_min(A)>=1. With approximate base solve, include the base residual: `||position residual||_2 <= R (||g||_2+||r_A||_2)` up to coordinate rounding. The old quantity R||g|| is a conservative **equation-residual** bound, not the true position error relative to a nonlinear root. F10J does not reinterpret the failed 1.850e-14 m bound as actual position error, and does not change the stopping rule. Independent root comparison supplies separate observed error evidence on this fixed state; no global a priori error theorem is claimed.

|Original linear solve / outer|Recursive norm|True norm|Original target|Decision|
|---|---:|---:|---:|---|
{linear_lines}

Base recursive/true residuals are {audit["linear"][0]["recursive_final"]:.9e} / {audit["linear"][0]["true_final"]:.9e}. They closely agree: this is an insufficient bounded iteration budget, not a catastrophic recursive-residual drift. Its independent position-equivalent linear error is only {audit["linear"][0]["max_position_equivalent_error_m"]:.6e} m, which also shows why a residual norm cannot be equated to solution error.

## Armijo diagnosis and minimal correction

The inherited rule is `E(q+t d) <= E(q)+1e-4 t g.d`, at most 20 halvings, at most 24 accepted updates. Final original g.d={last["direction_derivative"]:.9e}<0. All 20 failed candidates retain the same 38 active contacts. For the full candidate, direct float64 total-energy subtraction gives **+{last["original_energy_change"]:.9e}**; independent extended-precision stable difference gives **{last["stable_energy_change"]:.9e}**; Armijo requires at most {last["armijo_change_rhs"]:.9e}. The rejected full step has actual equation residual {last["equation_residual_m"]:.9e} m and original gradient bound {last["gradient_bound_m"]:.9e} m. Thus the stopping rule itself need not be loosened: the original code never accepted the available sufficient step.

`f10j_gpu.py` subclasses the entire original `FixedStateB2.solve`; it does not duplicate or replace the F10I solver. It reuses resistance, contact forces, normals, block preconditioning, PCG and original Armijo control flow. The bounded PCG budget becomes 256; true and recursive residual checks both gate admission. Conditional execution actually needs 132 base iterations and 93 for the final direction. No residual tolerance is loosened. A deliberately one-iteration base solve is rejected with TRUE_OR_RECURSIVE_LINEAR_RESIDUAL.

The merit is translated so current energy is zero and trial energy is the algebraically identical increment:

`Delta E = q^T A s + 0.5 s^T A s + alpha/2 sum (delta_new-delta_old)(delta_new+delta_old)`.

Distances use the rationalized increment `Delta d=(2 dx.ds+ds.ds)/(|dx+ds|+|dx|)`. For contacts staying active, `Delta delta=-Delta d`. This avoids subtracting two nearly equal total energies and two nearly equal distances. The inherited acceptance inequality then compares Delta E directly with 1e-4 t g.d. All {len(guards["merit_rows"])} original candidates agree with the independent long-double acceptance calculation in the actual GPU kernel check. Maximum absolute energy-increment discrepancy across the different-sized increments is {guards["maximum_merit_absolute_error"]:.3e}; no acceptance decision differs.

## Numerical evidence and scope

Original 18 small fixtures: 12 two-bubble and 6 four-bubble tests pass, including two expected overlap-budget rejections. The real 32-bubble case passes, as do nonfinite input, coincident-center, severe-overlap and incomplete-PCG rejection checks. Final conservative bound {new["conservative_equation_bound_m"]:.9e} m; base true residual {new["linear_audit"][0][1]:.9e}; overlap {m["max_overlap_m"] * 1e9:.9f} nm; maximum one-update correction {m["max_single_iteration_m"] * 1e9:.9f} nm. Original budgets 17.8/178 nm are unchanged. Contact/edge sets match CPU; reconstructed force mismatch {m["force_actual_geometry_error_N"]:.6e} N; net reaction {m["net_force_N"]:.6e} N. Nonnegative pair forces and zero force for separated pairs follow the unchanged geometric unilateral kernel and are audited against returned geometry.

GPU vs original CPU (atol=1e-14): position {m["max_position_vs_cpu_m"]:.9e} m; force {m["max_force_vs_cpu_N"]:.9e} N. GPU vs independent force-loop/Powell reference: position {m["max_position_vs_independent_m"]:.9e} m; force {m["max_force_vs_independent_N"]:.9e} N. Actual position equation residual {m["max_equation_residual_m"]:.9e} m. Both references solve this same exploratory discrete model; neither is experimentally calibrated truth.

## Synchronization semantics

Old intermediate reads: initial outer checks (5), after-direction checks (5), candidates (24), final check (1): 35. Corrected: 6 outer checks + 5 after-direction + 5 trials + final =17. Optimized reuses the **accepted** candidate residual for five immediately following outer checks:12. Rejected trials are never reused. Positions are bitwise identical across all 19 fixtures. Final host reads remain 12 old /14 new, including full vectors/operator evidence and both PCG checks; explicit prep/final synchronizations also remain.

Warp 1.18's actual `_run_capturable_loop` supports `capture_while`. Conditional bounded PCG graphs replace fixed 256-iteration captured loops without changing their criterion. Zero-RHS fixtures do not move; converged states are not subsequently changed. The original fixed-loop masked CG kernels also froze x after convergence, but still executed unused work. No new graph framework, GPU-resident nonlinear loop, batch integration or production default was added. Failed PCG is rejected at final admission; this fixed-state module is not a streaming accepted-step integrator.

Read-only instrumentation and ordinary original execution give exactly identical positions, also identical to archived F10I v2. Raw traces: `audit_v1/AUDIT.json`; independent audit: `analysis/ORIGINAL_AUDIT.json`; numerical runs: `corrected_v1`, `optimized_v2`; negative and merit checks: `guards_v1`.

Primary implementation references: [Warp 1.18 CG source](https://github.com/NVIDIA/warp/blob/v1.18.0/warp/_src/optim/linear.py), [SciPy CG stopping rule](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cg.html). Actual installed Warp source and SHA are archived in `audit_v1/WARP_CAPABILITIES.json`; project equations and independent checks remain the primary evidence.
"""
    (docs / "FAST_PZT_F10J_GPU_SOLVER_AUDIT.md").write_text(audit_text)
    costrows = "\n".join(
        f"|{n}|{medians[n]['median_wall_ms']:.6f}|{medians[n]['min_wall_ms']:.6f}–{medians[n]['max_wall_ms']:.6f}|{medians[n]['median_preparation_solve_ms']:.6f}|{medians[n]['median_quality_ms']:.6f}|{medians[n]['median_output_ms']:.6f}|"
        for n in names
    )
    profilerows = "\n".join(
        f"|{n}|{sum(v['graph_cuda_span_s'] for v in perf['diagnostics'][n]['pcg_graphs']) * 1e3:.6f}|{sum(v['wall_s'] for v in perf['diagnostics'][n]['evaluations'] if v['trial']) * 1e3:.6f}|{perf['diagnostics'][n]['result']['timing']['scalar_sync_read_wall_s'] * 1e3:.6f}|"
        for n in names[:3]
    )
    (docs / "FAST_PZT_F10J_PERFORMANCE.md").write_text(f"""# F10J fixed-state performance

**Stop GPU B2 production integration for this B=1, N=32 workload.** Numerical admission passed first; three matched measurements show optimized GPU {medians["optimized"]["median_wall_ms"]:.3f} ms versus CPU {medians["cpu"]["median_wall_ms"]:.3f} ms. GPU cost is {performance["new_gpu_over_cpu_cost"]:.3f} times CPU. This is a measured lack of advantage, not an extrapolated long-trajectory speed ratio.

Hardware: RTX 4090, Ryzen 9 7950X3D; Python 3.12.3, Warp 1.18.0+cu12, NumPy 2.2.6, SciPy 1.15.3. Existing isolated environment and driver unchanged. OPENBLAS/OMP/MKL/NUMEXPR threads=1 via the inherited wrapper. GPU ownership checked before each job; process monitoring archived. No unrelated GPU job was interrupted.

Same saved x/free, 10 us, R, viscosity, kc, discrete A, geometric forces, original quality budgets. CPU uses unchanged `geometric_contact` source with atol=1e-15, which additionally satisfies the original GPU 1e-14 gradient-bound target (verified), not merely its looser original CPU residual stop. Both methods meet the same independent position/force/geometry gates. CPU matrix assembly is **inside** timing. GPU coefficient/normal preparation, validation and complete final vector/operator/residual readback are inside. Same independent quality checks and same JSON byte payload fields are added to both. Disk evidence writing is outside both per-solve timers. Reference state loading is outside both. No matched output is skipped.

One warmup per method, then exactly three rotated-order pairs. First compile and reusable graph construction excluded from resident cost and recorded separately. One additional instrumented profile per GPU version is diagnostic only. Original F10I still fails and is not a valid production baseline. Historical 35.457 ms is not reused as a new result.

|Method|Full median ms|Range ms|Prep+solve median ms|Independent quality ms|Matched output ms|
|---|---:|---:|---:|---:|---:|
{costrows}

Independent component medians do not necessarily sum to the median whole measurement. All raw paired timing objects include CPU assembly/contact and GPU buffer preparation/nonlinear-loop/scalar waits/result readback. Construction records include buffer allocation, inherited graph build and replacement graph build; those inherited graphs remain part of this minimal subclass prototype's startup cost.

|Diagnostic profile|PCG graph device spans ms|Armijo evaluation wall ms|Scalar wait/read wall ms|
|---|---:|---:|---:|
{profilerows}

**These categories overlap.** In particular scalar read wall time includes waiting for earlier PCG work. It is not all PCIe copying or CPU overhead that can disappear by deleting a read. The solve-wide CUDA event span contains host decision gaps and is not pure kernel time. PCG-only events enclose graph execution without host decision gaps, but still include device graph/kernel scheduling overhead. The measured PCG work dominates this very small 96-DOF GPU solve; reducing 35 to 12 intermediate reads is insufficient.

New processes (one each, existing compiler cache, new CUDA context/buffers/graphs): {", ".join(r["method"] + " " + format(r["process_wall_s"], ".6f") + " s" for r in cold["rows"])}. Includes common harness imports, input loading, one solve/check/output, artifact write and shutdown. CPU imports the common Warp-containing harness but initializes no CUDA context. These are diagnostic startup costs, not a cold compiler benchmark or repeated throughput statistic.

All six registered remote wrappers total {sum(process_cost.values()):.6f} s; this is process wall time, not all development time or pure GPU occupancy. Independent child costs are already included in process_v1 and are not added again. No B64, 10 ms trajectory, AI dataset, pressure scan, new CFD or new acoustic solve was run.

Keep opt-in code for reproducibility and the CPU B2 as a finite-state numerical research reference. Do not continue production optimization merely because the old GPU failed or because scalar-read counts improved. A future architecture change would require a separate task, new admission and demonstrated workload need.
""")
    screen = js(BASE / "analysis/TOPOLOGY_SCREEN.json")
    candidate = js(BASE / "analysis/TOPOLOGY_CANDIDATE.json")
    dev = f"""# F10J — GPU B2 Solver Convergence and Execution Efficiency

The original GPU failure is explained and corrected without changing physics or tolerances. Corrected and optimized fixed-state solvers pass the registered numerical gates. **The optimized GPU is still slower than the matched CPU reference, so production integration is stopped.**

Baseline complete F10I commit: `{BASELINE}`; old GPU v2 execution `458dabada891b218d8dc91b256a6de6c3d7d496d`. Independent branch `dev/fast-pzt-f10j-gpu-b2-convergence-20261012`; WSL project `{ROOT}`. Remote `/workspace/mevo_acoustic_microbubbles/f10j`; immutable `.execution/<commit>` snapshots, unique Supervisor labels. F10I's {len(manifest)} manifest entries and all 13 original F5 input files are unchanged and hash-verified. Old compatibility symlinks and unrelated projects were not altered.

## Deliverables and actual execution

`src/mevo_mb/fast_sim/f10j_gpu.py` is an opt-in subclass, not a copy of the solver or production engine integration. `scripts/fast_pzt_f10j/` contains audit, numerical validation, fair timing, offline topology, figures and delivery scripts. Configurations separately preregister diagnosis, correction, optimization and timings. Complete results, negative cases, raw PCG/Armijo traces, logs, input/source hashes and actual execution IDs are in `results/fast_pzt/f10j/`. The final postprocessing commit is not substituted for a scientific execution commit.

"""
    for label, commit in gate["execution_commits"].items():
        dev += f"- `{label}`: `{commit}`.\n"
    dev += f"""
## Layered decisions

|Layer|Decision|
|---|---|
|Mathematical stopping, this fixed state|PASS; original 1e-14 m bound and original true/recursive PCG targets|
|GPU implementation, tested fixtures|PASS; 16 legal small states, 2 expected budget rejections, 1 real32; optimization bitwise equivalent|
|CPU/GPU performance|NO B=1 ADVANTAGE; GPU {medians["optimized"]["median_wall_ms"]:.3f} / CPU {medians["cpu"]["median_wall_ms"]:.3f} ms|
|Contact / near-field physical calibration|PHYSICS_UNCALIBRATED|
|Group position and shape|Prior F10H/I evidence retained, not newly expanded|
|Strict edge timing|Original 12/29 over 10 us, max130 us FAIL retained|
|Persistent split/merge prediction|PERSISTENT TOPOLOGY VALIDATION NOT RUN|
|AI data|training_ready=false; no production integration|

## Historical topology screening

All analysis uses existing trajectories, fixed IDs, original 2.1R graph and 2 ms dwell. {screen["reconstruction_status_counts"]["ACCEPTED_LEDGER_RECONSTRUCTION_PASS"]} historical ledgers reconstruct consistently. Six failed/terminated files have incomplete reconstruction coverage and are explicitly excluded; other skipped paths include empty schemas, multi-batch data and replicated throughput workloads. Skipped or incomplete data are not labelled as having no events. F10B has no additional trajectory.h5 in this archive; its references are inherited data. Full path/hash, model, accepted-frame scope, topology counts and exclusions are in TOPOLOGY_SCREEN.json. Histories are correlated runs, not independent samples.

F5 main contains 31 persistent aggregation merges before release, then 11 persistent splits and 2 merges after switch-off; one late split remains right-censored. F7 mature histories show no persistent events in examined records; F8 contains historical flow-release changes. F10D–H preparation variants have model-dependent merging. All three actual geometric B2 histories (F10H10/5 us, F10I tightened2ms) remain one connected group and provide no positive persistent-event accuracy test. Frozen-normal F10G B2, PBD M0 and geometric B2 are explicitly different models.

Candidate source: `{candidate["source"]}`, stored frame835 at417.5ms, all32 IDs active, finite positions, maximum overlap {candidate["max_overlap_m"] * 1e9:.6f}nm. Original M0 first persistent split at418.54ms. Proposed future examination window417.5–421.54ms. PZT/secondary off, original CFD on, original triangular wall; exact initial arrays and IDs in HISTORICAL_F5_CANDIDATE.npz, controls/inputs/hashes in TOPOLOGY_CANDIDATE.json. Minimum stored wall gap is approximately zero, so near-wall coupling and valid B2 restart must be checked first. This is a historical geometry candidate, **not an admitted B2 validation case or B2 truth**. No event is manufactured and no new topology trajectory was run.

## Reproduction and tests

Use the existing environment; remote jobs run committed snapshots through `python.sh` (BLAS threads1). Executed modules: `audit --label audit_v1`, `validate --label corrected_v1`, `validate --label optimized_v2 --optimized`, `performance --label performance_v1`, `cold_process --label process_v1`, `guards --label guards_v1`. Each label is immutable; choose a new label for a separately authorized rerun. Numerical tests run before optimization/performance. Old failure reproduction and GPU merit-kernel vs independent extended-precision tests are actual computations. No pytest that merely mirrors code is substituted for scientific checks.

Offline modules: `analyze_audit`, `topology_screen`, `figures`, `delivery`, `progress_document`, `verify`. Do not rerun scientific histories to regenerate plots. See the dedicated [solver audit](FAST_PZT_F10J_GPU_SOLVER_AUDIT.md) and [performance report](FAST_PZT_F10J_PERFORMANCE.md).

## Required Chinese answers

1. 原GPU卡在基础PCG、第一次方向PCG各自残差，以及最终原保守梯度界；末次Armijo把合格下降步拒绝了。
2. 递推与真实残差基本一致。基础为1.25499e−12／1.25371e−12，目标5.18976e−13；没有发现灾难性递推漂移。
3. 1.850e−14 m不是实际位置误差。旧程序停在尚未接受最后一个下降步的状态；不需要放宽1e−14门槛。
4. 两个接近的总能量计算出现舍入；完整候选被算成增加1.02349e−16，稳定差值实际减少约1.32813e−17，20个候选均被原比较拒绝。
5. 发现数值实现问题：基础及首方向未收敛仍继续、相近总能量比较不可靠。原物理方程及停止目标没有改错或放宽。
6. 修正后32泡通过全部登记的数学、参考和预算检查。保守界{new["conservative_equation_bound_m"]:.6e}m；五次全步更新。
7. 相对原CPU位置／力差{m["max_position_vs_cpu_m"]:.6e}m／{m["max_force_vs_cpu_N"]:.6e}N；相对独立参考{m["max_position_vs_independent_m"]:.6e}m／{m["max_force_vs_independent_N"]:.6e}N；真实残差{m["max_equation_residual_m"]:.6e}m。
8. 中途标量读回35→17→12；包含最终向量证据为47→31→26次数组读回。没有声称消除所有同步。
9. 本轮三组完整成本中位数：旧失败GPU{medians["original"]["median_wall_ms"]:.3f}ms、修正版{medians["corrected"]["median_wall_ms"]:.3f}ms、优化版{medians["optimized"]["median_wall_ms"]:.3f}ms。构图和独立进程另报，未冒用历史35.457ms。
10. 没有公平效率优势。CPU{medians["cpu"]["median_wall_ms"]:.3f}ms，优化GPU成本是其{performance["new_gpu_over_cpu_cost"]:.3f}倍。
11. 本轮停止GPU B2生产集成，保留独立模块和CPU有限状态研究参考，不继续昂贵长轨迹或大规模重写。
12. 找到F5关声后417.5ms候选，但源自M0且贴壁；需要另行验证B2合法重启，不能称已有B2持续分裂真值。
13. 原团级和逐泡位置／形状指标的短窗数值证据保留；严格连边时间仍失败，物理校准、持续拓扑预测与AI准入未完成。
14. corrected_v1与optimized_v2获固定状态数值准入；原F10I仍失败；负例正确拒绝；性能和生产准入未通过。
15. 四组实际英文PNG/PDF位于 `{BASE / "figures"}`；来源与人工阅图记录同目录。
16. 科研进度文档由独立脚本最小补充，备份为 `/home/lzy/projects/microbubble_aggregation.md.bak_F10J_20261012`；before/after、diff、SHA、标题结构检查保存在progress_document/。
17. F10K优先用CPU数值参考准备可复现的短窗持续拓扑与近壁耦合验证；先做重启和预算准入，不把M0事件时刻作为B2拟合目标，也不继续追逐130us瞬时边时差或未经验证的GPU生产吞吐。
"""
    (docs / "FAST_PZT_F10J_DEVELOPMENT.md").write_text(dev)
    print("Reports and layered FINAL_GATE generated; production integration stopped.")


if __name__ == "__main__":
    main()
