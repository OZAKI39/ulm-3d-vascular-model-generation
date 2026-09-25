"""Reference audit only. No hydrodynamic kernel, GCB closure, or simulation.

Historical files are read only. Public paper transcription candidates are
explicitly unverified; evaluating their printed expressions is a diagnostic,
not selecting them for the requested kernel.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    table = Path(args.table)
    expected = "71e62ef6bc22edc66ae0c4d1dcc224378c2f6c5d3848578f3fbedd2cd9a43f98"
    assert hashlib.sha256(table.read_bytes()).hexdigest() == expected, "STOP_STATE_DRIFT"
    with h5py.File(table, "r") as f:
        attrs = {k: (v.item() if hasattr(v, "item") else v) for k, v in f.attrs.items()}
        datasets = {k: {"shape": list(f[k].shape), "dtype": str(f[k].dtype),
                        "attributes": dict(f[k].attrs)} for k in f}
        eps = f["epsilon"][:]
        excess = f["R_wall_excess_scaled"][:]
        total = f["R_total_scaled"][:]
        total_si = f["R_total_RMBW_SI"][:]
        flags = [json.loads(x) for x in f["reference_mode_flags"][:]]
    assert np.isfinite(excess).all() and np.isfinite(total).all()
    assert np.all(np.diff(eps) > 0)
    a, mu = attrs["reference_radius_m"], attrs["reference_mu_Pa_s"]
    scale = np.sqrt([6 * np.pi * mu * a] * 3 + [8 * np.pi * mu * a**3] * 3)
    identity_error = float(np.max(np.abs(total - excess - np.eye(6))))
    scaling_error = float(np.max(np.abs(total_si / np.outer(scale, scale) - total)) / np.max(np.abs(total)))
    assert identity_error <= 1e-12, "STOP_REFERENCE_AMBIGUITY"
    assert scaling_error <= 1e-12, "STOP_UNIT_AMBIGUITY"
    metadata = {
        "table_sha256": expected, "size_bytes": table.stat().st_size,
        "attributes": attrs, "datasets": datasets,
        "epsilon_domain": [float(eps[0]), float(eps[-1])], "grid_points": len(eps),
        "stored_runtime_matrix": "R_wall_excess_scaled = R_total_scaled - I",
        "total_minus_excess_minus_identity_maxabs": identity_error,
        "SI_reconstruction_relative_error": scaling_error,
        "radius_viscosity": "Dimensionless epsilon table; reference a,mu are metadata. S=diag(sqrt(6*pi*mu*a) x3,sqrt(8*pi*mu*a^3) x3); R_excess_SI=S E S.",
        "DOF_order": ["V_t1", "V_t2", "V_n", "Omega_t1", "Omega_t2", "Omega_n"],
        "load_order": ["F_t1", "F_t2", "F_n", "T_t1", "T_t2", "T_n"],
        "load_sign": "External applied load; torque about sphere center",
        "frame": "Right-handed Q=[t1 t2 n], n wall->center; T=diag(Q,Q); global=T local T^T",
        "historic_interpolation": "Full 6x6 componentwise linear interpolation in log(epsilon)",
        "historic_exterior_behavior": "epsilon<=0.001 clamp; epsilon>20 zero excess. NOT adopted for Phase 1.",
        "required_kernel_exterior_behavior": "h<=0 INVALID_GAP; no silent table extrapolation",
        "limitations": json.loads(attrs["limitations"]),
        "classification": "REFERENCE_WITH_LIMITATIONS",
        "table_min_symmetric_total_eigenvalue": float(np.linalg.eigvalsh((total + total.transpose(0, 2, 1))/2).min()),
        "eigenvalue_scope": "Input table diagnostic only; not a C++ kernel SPD test",
        "mode_flags_examples": {str(float(eps[j])): flags[j] for j in [0, len(eps)//2, len(eps)-1]},
    }
    write_json(root / "reference/RMBW_TABLE_METADATA.json", metadata)
    with (root / "raw/resistance_reference.csv").open("w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["epsilon", "normal_TT_over_6pi_mu_a", "tangential_TT_over_6pi_mu_a",
                         "parallel_RR_over_8pi_mu_a3", "normal_RR_over_8pi_mu_a3",
                         "TR_x_Omega_y_over_6pi_mu_a2"])
        for e, r in zip(eps, total):
            writer.writerow([e, r[2, 2], r[0, 0], r[3, 3], r[5, 5], r[0, 4]*math.sqrt(4/3)])

    # Small public-transcription candidates from GCB p657, not authenticated data.
    candidates = [("Table2", math.cosh(.08)-1, .45291, .48300),
                  ("Table2", math.cosh(.1)-1, .47861, .50818),
                  ("Table2", math.cosh(.3)-1, .65375, .67462),
                  ("Table3", .001, .3966, .4268),
                  ("Table3", .0001, .3185, .3468)]
    rows = []
    for origin, e, u, w in candidates:
        den = .6376 - .200*math.log(e)
        fu, fw = .7431/den, .8436/den
        rows.append({"table": origin, "epsilon": e, "transcribed_FU": u,
                     "transcribed_FOmega": w, "printed_expression_FU": fu,
                     "printed_expression_FOmega": fw,
                     "FU_relative_difference_to_transcribed_value": abs(fu/u-1),
                     "FOmega_relative_difference_to_transcribed_value": abs(fw/w-1),
                     "qualification": "UNVERIFIED_TRANSCRIPTION_DIAGNOSTIC",
                     "enters_physics_pass": False})
    with (root / "raw/gcb_transcription_differences.csv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_json(root / "reference/GCB_REFERENCE_AUDIT.json", {
        "requested_primary_doi": "10.1016/0009-2509(67)80048-4",
        "publisher_url": "https://www.sciencedirect.com/science/article/pii/0009250967800484",
        "primary_page_transcription": "https://www.scribd.com/document/684300240/Goldman1967-b",
        "locator": "Printed p657, Eqs (4.11),(4.12), Tables 2 and 3",
        "primary_page_image_verified": False, "trusted_local_numerical_reference": False,
        "printed_candidates": {"FU": "0.7431/(0.6376-0.200*ln(epsilon))",
                               "FOmega": "0.8436/(0.6376-0.200*ln(epsilon))"},
        "formula_class_if_resolved": "NEAR_WALL_ASYMPTOTIC_ONLY",
        "normalization": "Paper center height h -> project d; paper surface gap delta -> project h. FU=U/(d*S), FOmega=Omega/(S/2).",
        "asymptotic_domain": "epsilon << 1; no certified finite interval established",
        "later_paper_doi": "10.1093/qjmam/56.3.381",
        "later_paper_transcription": "https://doczz.net/doc/7418969/creeping-flow-around-a-sphere-in-a-shear-flow-close-to-a-...",
        "later_paper_locator": "Section 6, printed p404, Eqs (6.5),(6.6) and discussion following (6.8)",
        "later_paper_observation": "Discusses poor accuracy of the printed translation expression and a possible misprint. Its own approximants were not adopted.",
        "later_paper_corrigendum_doi": "10.1093/qjmam/hbs012",
        "corrigendum_checked": False,
        "modern_reproduction": "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023WR035456",
        "modern_locator": "Table 1 F3/F6 uses these expressions for H<0.05; does not resolve their independent accuracy or source discrepancy",
        "finite_gap_interpolation_selected": False,
        "candidate_table_points": rows,
        "decision": "STOP_GCB_REFERENCE_INSUFFICIENT",
        "reason": "No authenticated primary-page numerical dataset; printed-formula versus transcribed-table discrepancy unresolved. Cannot select or repair a continuous closure without guessing.",
        "qualification": "Evidence of an unresolved reference discrepancy, not proof that the primary paper itself is erroneous.",
    })

    source_base = "/home/lzy/projects/compre_output/wall_hydrodynamics_reference_audit/20260916_115322/src/"
    capability = {}
    specs = {
        "normal_translation": ("classical_wall_reference.py::normal_brenner", "0.001..20 tabulated; prior 14-gap Brenner comparisons", "SUPPORTED_HISTORICALLY; not requalified"),
        "tangential_translation": ("compare_oneill_table.py; classical_wall_reference.py", "0.001..20 table; O'Neill discrete alpha points; asymptotic epsilon<<1", "SUPPORTED_WITH_DOMAIN_LIMITS"),
        "parallel_rotation": ("classical_wall_reference.py::published_wall_asymptotics", "0.001..20 table; asymptotic epsilon<<1", "PARTIAL; finite-gap bound and 0.3817/0.3709 constant unresolved"),
        "translation_rotation": ("compare_oneill_table.py; classical_wall_reference.py", "0.001..20 table; O'Neill discrete points; near-wall author-fitted linear term", "PARTIAL; continuous-gap accuracy unverified"),
        "normal_axis_rotation": ("classical_wall_reference.py::published_wall_asymptotics", "0.001..20 table; contact-limit support only", "PARTIAL; finite-gap accuracy unverified"),
    }
    for name, (src, domain, status) in specs.items():
        capability[name] = {"available": True, "source": source_base + src, "gap_domain": domain,
                            "availability_meaning": "Historical coefficient/reference data exist; not all-gap exact accuracy",
                            "qualification": status}
    for mode in ["translation", "rotation"]:
        capability["gcb_free_shear_"+mode] = {
            "available": False, "source": "No closure in the three permitted historical scripts; see reference/GCB_REFERENCE_AUDIT.json",
            "local_implementation_available": False,
            "public_asymptotic_expression_found": True,
            "qualified_continuous_reference_available": False,
            "gap_domain": "NONE_CERTIFIED",
            "qualification": "STOP_GCB_REFERENCE_INSUFFICIENT"}
    capability["resistance_table"] = metadata
    capability["historical_tolerances"] = {
        "wall_table_interpolation_action_relative": .001,
        "wall_table_reciprocity_relative": 1e-10,
        "historical_close_agreement_label": "<0.02; descriptive, NOT physical acceptance",
        "historical_physical_accuracy_gate": "No universal 1% or agreement-based pass",
        "GCB_physical_tolerance": None,
        "recommendation": "Establish a reference uncertainty budget and approve translation/rotation error limits before physical certification; no numerical tolerance invented here."}
    write_json(root / "REFERENCE_CAPABILITY.json", capability)
    tests = ["matrix_symmetry", "positive_dissipation_SPD", "SI_scaling", "frame_rotation_invariance",
             "tangent_basis_sign_flip", "zero_shear", "shear_sign_reversal", "shear_magnitude_linearity",
             "no_fake_normal_shear_forcing", "isolated_sphere_GCB_linear_solve_recovery"]
    write_json(root / "validation/PHASE1_VALIDATION.json", {
        "status": "BLOCKED", "stop_code": "STOP_GCB_REFERENCE_INSUFFICIENT",
        "algebra_pass": False, "algebra_status": "NOT_RUN_REFERENCE_STOP",
        "physics_reference_status": "BLOCKED", "certified_epsilon_domain": [],
        "standalone_kernel_created": False, "independent_python_kernel_validator_created": False,
        "tests": {f"{chr(65+i)}_{name}": "NOT_RUN_REFERENCE_STOP" for i, name in enumerate(tests)},
        "input_table_checks": {"total_equals_bulk_plus_excess": identity_error,
                                "SI_reconstruction_relative_error": scaling_error},
        "input_checks_are_kernel_tests": False,
        "gcb_reference_diagnostics": rows,
        "shear_rhs_components_plot": "NOT_GENERATED_NO_QUALIFIED_CLOSURE",
        "invalid_gap_NaN_Inf_tests": "NOT_RUN_NO_KERNEL",
        "no_scientific_pass_claim": True})
    print(json.dumps({"status": "BLOCKED", "stop_code": "STOP_GCB_REFERENCE_INSUFFICIENT",
                      "table_identity_error": identity_error, "SI_scaling_error": scaling_error,
                      "epsilon_0p001_diagnostic": rows[3]}, indent=2))


if __name__ == "__main__":
    main()
