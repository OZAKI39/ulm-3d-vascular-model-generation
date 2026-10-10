"""Evidence-linked Stage 4 report, including the Stage I.5 performance audit."""
from pathlib import Path
from html import escape
from scripts.stage4.common import ROOT,new_run,require,write,sha,finish
from scripts.stage4.postprocess import read,stage,secondary_comparison,performance,paraview,saved_slip_reynolds
from scripts.stage4.figures import make_figures
from scripts.stage4i5.performance_figures import make_performance_figures


def number(value,digits=4):
    return "undefined" if value is None else f"{value:.{digits}g}"


def interval(metric):
    ci = metric["bootstrap_95CI"]
    return number(metric["mean"])+(" ["+number(ci[0])+", "+number(ci[1])+"]" if ci is not None else " [undefined]")


def table(headers,rows):
    return "| "+" | ".join(headers)+" |\n|"+"|".join(["---"]*len(headers))+"|\n"+"\n".join("| "+" | ".join(map(str,row))+" |" for row in rows)


def link(path,title=None):
    return f"[{title or path}](../{path})"


def main():
    stages = {name:stage(name) for name in ["audit","audit/final","geometry","synthetic_flow","simvascular/steady",
        "simvascular/import","acoustic","phase_screening","primary_only","statistics","statistics/paired_runs",
        "primary_secondary","convergence","history_branch","performance_refactor/backend_decision"]}
    require("audit/final")
    for name in stages:
        if name != "history_branch":
            require(name)
    history = stages["history_branch"]
    if history["status"] not in ["PASS","PASS_WITH_SCOPE","WARNING"]:
        raise RuntimeError("Incomplete or failed history reference must be resolved/documented before final report")
    out = new_run("report",dict(stage_runs={name:s["run"] for name,s in stages.items()}),["audit/final"])
    secondary = secondary_comparison(out)
    perf = performance(out)
    collections = paraview(out)
    figures = make_figures(out,secondary)
    decision=stages["performance_refactor/backend_decision"]
    figures.extend(make_performance_figures(out,decision))
    stats,conv = stages["statistics"],stages["convergence"]
    effects = stats["paired_effects"]
    mirror = (stats["cases"]["B2"]["f_A"]["mean"] is not None and stats["cases"]["B3"]["f_A"]["mean"] is not None
              and stats["cases"]["B2"]["f_A"]["mean"]>.5 and stats["cases"]["B3"]["f_A"]["mean"]<.5
              and effects["B2_minus_B3_f_A"]["bootstrap_95CI"] is not None and effects["B2_minus_B3_f_A"]["bootstrap_95CI"][0]>0)
    lower = conv["paired_effect_lower_95CI"]
    uncertainties = conv["observed_fraction_uncertainties"]
    robust = {name:lower is not None and error is not None and lower>error for name,error in uncertainties.items()}
    robust["timestep"] &= conv["timestep_outcomes_robust"]
    grid_outcome=conv["acoustic_grid"]["branch_comparison"]
    robust["grid"] &= (grid_outcome["assignment_agreement"]>=.95 and
        max(abs(grid_outcome["A_over_initial_difference"]),abs(grid_outcome["B_over_initial_difference"]))<=.05)
    items = {}
    def item(name,status,evidence,scope):
        paths = [str(Path(p)) for p in evidence]
        items[name] = dict(status=status,evidence=[dict(path=p,sha256=sha(ROOT/p)) for p in paths],scope=scope)
    for name,key in [("geometry","geometry"),("synthetic_flow","synthetic_flow"),("simvascular_flow","simvascular/steady"),
                     ("flow_import","simvascular/import"),("acoustic_basis","acoustic")]:
        s = stages[key]
        item(name,s["status"],[str(Path(s["run"])/"SUMMARY.json")],s.get("scope","Recorded checks"))
    item("branch_classifier","PASS_WITH_SCOPE",[str(Path(stages["audit/final"]["run"])/"SUMMARY.json"),
         str(Path(stages["synthetic_flow"]["run"])/"SUMMARY.json")],"First outward finite-plane crossing, immutable global IDs; geometric tolerance diagnostic reported separately")
    item("branch_control","PASS_WITH_SCOPE" if stats["statistical_directional_control_supported"] else "WARNING",
         [str(Path(stats["run"])/"SUMMARY.json")],"Paired-seed conditional exit fractions; raw delivery and unassigned outcomes remain distinct. Combined phase/amplitude control, not phase-only attribution")
    item("mirror_control","PASS_WITH_SCOPE" if mirror else "WARNING",[str(Path(stats["run"])/"SUMMARY.json")],
         "B2 prefers A, B3 prefers B and paired interval separates them; no exact particle-level symmetry required")
    item("secondary_effect","PASS_WITH_SCOPE",[str((out/"SECONDARY_COMPARISON.json").relative_to(ROOT))],secondary["scope"])
    for name,key in [("timestep_robustness","timestep"),("grid_robustness","grid")]:
        item(name,"PASS_WITH_SCOPE" if robust[key] else "WARNING",[str(Path(conv["run"])/"SUMMARY.json")],
             "One seed/control, two resolutions. Conditional-fraction difference is compared with paired-effect lower CI; per-ID agreement >=95% and initial-population arrival differences <=5 percentage points are additional declared diagnostic gates. Timestep early-cluster relative L2 gate is 10%. These are reporting policies, not proof of asymptotic convergence.")
    item("history_branch_robustness",history["status"],[str(Path(history["run"])/"SUMMARY.json")],history["scope"])
    item("performance","PASS_WITH_SCOPE",[str((out/"PERFORMANCE.json").relative_to(ROOT)),
         str(Path(decision["run"])/"SUMMARY.json"),"results/stage4/performance_refactor/gold_baseline/GOLD_BASELINE_MANIFEST.json"],
         "Immutable full gold plus matched short replay gates; algorithmic and geometry correctness precedes measured speed. The 2x target is not a scientific pass criterion.")
    model_evidence = []
    for row in stages['primary_only']['rows']+stages['primary_secondary']['rows']:
        summary = read(Path(row['run'])/'particles/SUMMARY.json')
        validity = summary['physics_validity']
        model_evidence.append(dict(case=row['case'],secondary=row['secondary'],run=row['run'],
            validity_status=validity['validity_status'],physics_status=validity['physics_status'],
            wall_hit_particle_count=summary['wall_hit_particle_count'],
            invalid_accepted_duration_s=validity['invalid_accepted_duration_s'],
            accepted_duration_s=validity['accepted_duration_s'],
            saved_slip_reynolds=saved_slip_reynolds(row),
            duration_semantics="Duration whose accepted endpoint has at least one validity reason; not a fraction of individual trajectories"))
    finite_re_secondary = any(r['secondary'] and r['saved_slip_reynolds']['measurements']['constrained']['maximum'] is not None
        and r['saved_slip_reynolds']['measurements']['constrained']['maximum']>=1 for r in model_evidence)
    if finite_re_secondary:
        items['secondary_effect']['status'] = 'WARNING'
        items['secondary_effect']['scope'] += ' Saved constrained-slip Re reaches/exceeds 1; the frozen Stokes approximation is outside its small-Re regime. Quantification is descriptive model extrapolation.'
    success = bool(stats["statistical_directional_control_supported"] and mirror and all(robust.values())
        and history["status"] in ["PASS","PASS_WITH_SCOPE"] and not finite_re_secondary)
    acceptance = dict(status="PASS_WITH_SCOPE" if success else "WARNING",
        stage_state="BIFURCATION_CONTROL_VALIDATED_WITH_SCOPE" if success else
            "BIFURCATION_CONTROL_DEMONSTRATED_WITH_QUANTITATIVE_WARNING" if stats["statistical_directional_control_supported"] else "BIFURCATION_CONTROL_NOT_VALIDATED",
        items=items,stage3_baseline_integrity=stages["audit/final"],same_initial_particles_verified=stats["paired_identity_verified"],
        production_history_free=not stages["audit/final"]["production_history_imports"],flowfield_interface="velocity(x, t)",
        patient_specific_started=False,statistics=stats["run"],numerical_convergence=conv["run"],
        model_validity=model_evidence,near_wall_hydrodynamics_validated=False,
        finite_Re_secondary_warning=finite_re_secondary,finite_Re_corrections_in_production=False,
        concentration_claim="Normalized finite-ensemble occupancy only; no clinical injection concentration",
        control_causality="Selected controls change both phase and amplitude; isolated phase causality was not tested",
        report_run=str(out.relative_to(ROOT)))
    write(out/"STAGE4_ACCEPTANCE.json",acceptance)
    write(ROOT/"STAGE4_ACCEPTANCE.json",acceptance)
    primary = {r["case"]:r for r in stages["primary_only"]["rows"]}
    cfd,imp,acoustic = stages["simvascular/steady"],stages["simvascular/import"],stages["acoustic"]
    sections = []
    def section(title,text):
        sections.append(f"## {len(sections)+1}. {title}\n\n{text}\n")
    transition_tables = []
    outcomes = ["BRANCH_A","BRANCH_B","UNASSIGNED","WALL_LOSS"]
    for control in ["B2","B3"]:
        counts = stats["pooled_transition_counts"]["B0->"+control]
        transition_tables.append(f"**B0 → {control} matched-particle transition counts, pooled across seeds**\n\n"+
            table(["B0 outcome / control outcome","A","B","U","W"],
                  [[label]+[counts[origin][destination] for destination in outcomes]
                   for label,origin in zip(["A","B","U","W"],outcomes,strict=True)]))
    section("Executive summary",f"Stage state: **`{acceptance['stage_state']}`**. All requested Stage 4A–K executions and analyses are recorded. "
        f"Paired primary-only directional control supported: **{stats['statistical_directional_control_supported']}**; "
        f"the paired conditional-exit-fraction effect exceeds the tested conditional-fraction differences: **{conv['effect_exceeds_tested_numerical_uncertainties']}**. "
        f"The additional per-particle and absolute-arrival gates give timestep robustness **{robust['timestep']}** and grid robustness **{robust['grid']}**. "
        f"The selected history-reference branch agreement is {history['branch_agreement_rate']:.1%}.\n\n"
        "The outcome is a 6 s finite initial-ensemble experiment, with 256 bubbles per seed, 10 paired seeds, "
        "2.5 MHz homogeneous-water acoustics and a steady SI SimVascular Y flow. Unassigned bubbles are observed without a classified exit by 6 s. "
        "Conditional branch enrichment does not imply increased absolute delivery to the target. Controls vary both phase and amplitude.\n\n"+
        "Near-wall physics remains outside the validated model domain. In the first primary-only seed, "
        +", ".join(f"{r['case']}: {r['wall_hit_particle_count']}/256 particles projected at a wall" for r in model_evidence if not r['secondary'])+
        ". Agreement between numerical resolutions or history references does not validate the omitted wall hydrodynamics.\n\n"+
        f"Secondary finite-Re warning: **{finite_re_secondary}**. A measured constrained-slip Re ≥ 1 is outside the small-Re Stokes assumption; completing an extrapolative calculation does not establish its physical validity.\n\n"+
        link("STAGE4_ACCEPTANCE.json","Machine-readable acceptance"))
    section("Stage3 baseline integrity",f"Frozen tag `stage3-minimal-validated` points to `{stages['audit']['baseline_commit']}`. "
        f"The final audit rechecked {stages['audit/final']['protected_files']} original files; changed: {stages['audit/final']['changed']}. "
        f"It also verified the hashes of {stages['audit/final']['baseline_metadata_files']} original acceptance/report files against the frozen commit. "
        f"It also verified that {stages['audit/final']['stage4_unchanged_gold_physics_source_count']} original production physics source files remain identical across {stages['audit/final']['stage4_completed_run_count']} completed Stage 4 particle runs. Added optional performance modules are audited by their versioned per-run source hashes and the I.5 gate. "
        "Original PHYSICS_ACCEPTANCE and accepted Stage 3 results retain their identities.\n\n"+
        "\n\n".join(f"**{c['name']}**: return code {c['returncode']}.\n```text\n{c['output'].strip()}\n```" for c in stages["audit/final"]["checks"])+"\n\n"+
        link(str(Path(stages['audit/final']['run'])/'SUMMARY.json'),"Final integrity and regression evidence")+"; "+
        link(str(Path(stages['audit/final']['run'])/'STAGE3_BASELINE_METADATA.json'),"Frozen acceptance/report hashes"))
    section("Bifurcation geometry","Config-driven smooth symmetric Y: parent 14 mm / radius 1.3 mm; daughters 12 mm / radius 1 mm, angles ±30°, smooth transition scale 0.2 mm. "
        "The positive-y triangulated surface half is reflected for exact surface symmetry. Meshing introduces a small asymmetric triangulation difference; the CFD split is checked independently. "
        f"Geometry hash `{stages['geometry']['geometry_hash']}`. Surface is closed and connected; native tetrahedral mesh has 26,712 nodes and 146,230 noninverted tetrahedra.\n\n"+
        link(stages['geometry']['surface'],"Labelled surface")+"; "+link("configs/stage4/bifurcation.yaml","Geometry/flow/initialization configuration"))
    section("Branch classification","Labels: PARENT, BRANCH_A, BRANCH_B, UNASSIGNED, WALL_LOSS. Classifier checks every accepted substep for the first outward crossing of a finite daughter plane 10 mm from the bifurcation origin. "
        "Time, position, velocity, transverse velocity and global particle ID are interpolated at that event. Assignment is permanent. Exited bubbles are removed from force/flow queries; HDF5 inactive positions remain frozen at the measured crossing, marked inactive, with inactive velocities/forces NaN. "
        "Contact or wall projection is not automatically WALL_LOSS. Saved visualization frames are every 20 ms; classifier events are not downsampled. Initial tolerance 1 nm; sensitivity uses 0.1 µm. "
        "Bifurcation-plane positions are separate first-crossing events, never guessed from terminal coordinates.")
    section("Synthetic-flow software check", "`SOFTWARE_TEST_ONLY`: smooth bounded divergence-free streamfunction transport checks routing, first crossing and paired bookkeeping. "
        "At the accepted 12 s software observation window, uniform and flux-weighted initializations both classify all 256 particles; uniform A/B=128/128, flux A/B=130/126. "
        "An earlier 3 s check with many unassigned particles remains archived; its observation limit was insufficient. Synthetic flow is excluded from all scientific branch-control claims.\n\n"+
        link(str(Path(stages['synthetic_flow']['run'])/'SUMMARY.json'),"Synthetic software validation"))
    section("SimVascular CFD","Existing SimVascular/TetGen and the existing GPU-enabled svMultiPhysics executable were reused; no CFD solver was implemented. "
        "Water density 1,000 kg/m³, viscosity 0.001 Pa·s; parabolic inlet mean speed 0.005 m/s; equal 0 Pa reference-pressure Neumann outlet conditions; no-slip wall. These are traction conditions, not an imposed uniform native pressure array. "
        "Five consecutive interval checks reach E_u<1e−5 and E_Q<1e−6 at step 120; native final output is step 121.\n\n"+
        table(["Quantity","Measured value"],[["Q inlet [m³/s]",number(cfd['metrics']['Q_in_m3_s'],10)],
            ["Q A [m³/s]",number(cfd['metrics']['outlet_flows_m3_s']['OUTLET_A'],10)],
            ["Q B [m³/s]",number(cfd['metrics']['outlet_flows_m3_s']['OUTLET_B'],10)],
            ["Relative mass imbalance",number(cfd['metrics']['epsilon_mass'])],
            ["A fraction of CFD flow",number(cfd['metrics']['outlet_fractions']['OUTLET_A'],10)],
            ["Peak speed [m/s]",number(cfd['metrics']['velocity_max_m_s'])],
            ["Wall speed maximum [m/s]",number(cfd['metrics']['wall_velocity_max_m_s'])]])+
        "\n\nCoordinates are metres, velocity m/s, time seconds, pressure Pa by generated-input/solver contract. This checks steady numerical coupling. Global conservation and GPU equivalence do not establish near-wall/junction local physical accuracy; CFD mesh independence was not tested.")
    section("Flow import validation",f"The existing `velocity(x, t)` interface is unchanged. Native array `{imp['velocity_array']}`, coordinate and velocity scales both 1. "
        "Default VTK probing could return a neighboring tetrahedron for a few near-face coordinates even at zero tolerance. The Stage 4 option `tetra_p1` retains a cached VTK locator but verifies barycentric membership and checks exact candidate tetrahedra before P1 interpolation. Outside points raise an error; no nearest-neighbor extrapolation. Default Stage 3 behavior remains unchanged.\n\n"+
        table(["Sample group","N","Max velocity error [m/s]","Relative L2"],[[key,row['count'],number(row['max_absolute_error_m_s']),number(row['relative_L2'])] for key,row in imp['comparisons'].items()])+
        f"\n\nSnapshot SHA-256 `{imp['snapshot_sha256']}`. Batch-256 interpolation median {imp['batch256_median_s']:.6g} s. "
        "These comparisons establish import fidelity to the native discrete solution, not physical accuracy.\n\n"+link(str(Path(imp['run'])/'SIMVASCULAR_FLOW_IMPORT_VALIDATION.json'),"Import validation"))
    controls = stages['phase_screening']['cases']
    section("PZT acoustic basis",f"New geometry-keyed basis SHA-256 `{acoustic['basis_sha256']}`; four 6 mm ideal acoustic discs around the bifurcation. "
        "k-Wave on RTX4090 uses 0.1 mm grid, 256×192×128 total grid; saved complex ROI covers distal parent, bifurcation and daughters. "
        f"Independent direct superposition errors: pressure {acoustic['pressure_superposition_relative_L2']:.4g}, gradient {acoustic['gradient_superposition_relative_L2']:.4g}. "
        "Phasor convention is Re[P exp(+iωt)]. A field-only coarse screen selected candidate 42 before formal particle trials; the eligible controls use all four elements and nonzero relative y-element phase. This is not global optimization.\n\n"+
        table(["Control","Amplitudes [Pa], ordered y+, y−, z+, z−","Phases [rad]"],[[label,", ".join(number(x,7) for x in controls[label]['amplitudes_Pa']),", ".join(number(x,5) for x in controls[label]['phases_rad'])] for label in controls])+
        "\n\nThe Y field includes pressures below the former 5 kPa lower table bound. A separate offline adaptive table extends the same frozen APECSS/Marmottant/RPAR model to 0.1–40,000 Pa on the narrow R=1.779–1.781 µm / f=2.499–2.501 MHz manifold; 128 independent holdouts pass. The original Stage 2.5 table is unchanged. The declared effective radial viscosity is 0.002 Pa·s; CFD and translation viscosity remain 0.001 Pa·s, preserving Stage 3's separation of these parameters. There is no runtime pressure clipping/extrapolation.")
    section("Primary-only branch steering","Formal seed 20261002 uses flux-weighted sampling at x=−2 mm within 90% of the parent radius, fixed IDs 0–255, R0=1.78 µm, minimum initial spacing 8 µm. Every control reuses exactly identical initial coordinates/radii/IDs. Acoustic-off B0 bypasses the positive-pressure table correctly.\n\n"+
        table(["Case","A","B","U","Wall loss","f_A among exits","Runtime [s]"],[[label,row['N_A'],row['N_B'],row['N_unassigned'],row['N_wall_loss'],number(row['f_A']),number(row['elapsed_s'],6)] for label,row in primary.items()])+
        "\n\nU is observation censoring, not an accepted exit or a numerical failure. Old failed and interrupted development attempts remain identifiable; they are not included in the formal completed-run table.")
    section("Mirror phase control",f"Mirror-control criterion supported: **{mirror}**. B3 swaps the two y-element complex controls of B2, retaining the z pair. "
        "The field mirror error is recorded in the screen evidence; finite random initial populations are shared across controls, not reflected individually. "
        f"Paired B2−B3 f_A mean and CI: {interval(effects['B2_minus_B3_f_A'])}. Mirror residual f_A(B2)+f_A(B3)−1: {interval(effects['mirror_fA_sum_minus_1'])}. "
        "This comparison concerns the combined phase/amplitude patterns; it does not isolate phase as the sole causal parameter.")
    section("Multi-seed statistics",f"{stats['seed_count']} independent paired seeds: {stats['seeds']}. For each seed, all four controls share initial arrays by exact byte hash. "
        "20,000 bootstrap samples resample seed blocks; particles within a run are not treated as independent seed replicates. No seed is dropped when a conditional fraction is undefined.\n\n"+
        table(["Case","Mean f_A [95% CI]","SD f_A","Mean A/N0 [95% CI]","Mean B/N0 [95% CI]","Mean U/N0 [95% CI]"],
              [[label,interval(s['f_A']),number(s['f_A']['std']),interval(s['A_fraction_of_initial']),interval(s['B_fraction_of_initial']),interval(s['unassigned_fraction'])] for label,s in stats['cases'].items()])+
        "\n\n"+table(["Paired effect","Mean [95% CI]"],[[key,interval(value)] for key,value in effects.items()])+
        "\n\n"+link(str(Path(stats['run'])/'BRANCH_SWITCH_TABLE.csv'),"Complete A/B/U/W transition counts")+"; "+link(str(Path(stats['run'])/'PARTICLE_OUTCOMES.csv'),"Crossing/residence/censoring/bifurcation outcomes")+
        ". Residence time starts at the distal-parent initialization plane (x=−2 mm), not at the upstream CFD inlet (x=−14 mm). Right-censored residence times are observation durations, not known final residence times. Proximal ROI normalized occupancy and target/non-target concentration ratios are in each run's metrics.json; zero denominators remain undefined.\n\n"+
        "\n\n".join(transition_tables)+"\n\nPooled transition counts describe actual per-ID switches and censoring. Inferential intervals above resample whole paired seeds, not individual particles.")
    secondary_outcomes = []
    for label,s in secondary['cases'].items():
        target = 'f_A' if label=='B2' else 'f_B'
        for mode,title in [('primary','Primary'),('secondary','Primary + secondary')]:
            c = s[mode+'_counts']
            secondary_outcomes.append([label,title,c['N_A'],c['N_B'],c['N_unassigned'],c['N_wall_loss'],
                c['N_initial'],number(c[target]),number(c['target_fraction_of_initial_population'])])
    section("Secondary interaction effect",table(["Case","Physics","A","B","U","W","N0","Conditional target f","Target / N0"],secondary_outcomes)+
        "\n\nU denotes no classified exit by the observation cutoff; W denotes wall loss. Conditional target f uses A+B as its denominator; Target/N0 uses the complete initial population.\n\n"+
        table(["Case","Δ target f (secondary − primary)","Δ target/N0 (secondary − primary)","Changed assignments","NN early slope primary [m/s]","NN early slope secondary [m/s]","Secondary near-contact maximum"],
        [[label,number(s['target_fraction_difference']),number(s['target_arrival_fraction_difference']),s['changed_assignments'],number(s['primary']['early_NN_slope_m_s']),number(s['secondary']['early_NN_slope_m_s']),number(s['secondary']['near_contact_fraction']['maximum'])] for label,s in secondary['cases'].items()])+
        "\n\n"+secondary['scope']+" Detailed mean NN, 0–0.1 s NN slope, cluster count, largest cluster fraction, particle-contact fraction, far/transition/near/contact pair fractions and A/B spatial cluster metrics are in "+link(str((out/'SECONDARY_COMPARISON.json').relative_to(ROOT)),"SECONDARY_COMPARISON.json")+". The finite surviving population changes with exit; clustering is geometric, not coalescence. No assumption that secondary interactions improve steering is made.")
    section("Numerical convergence",table(["Comparison","Assignment agreement","Δ f_A","Δ A/N0","Δ B/N0","Mean |Δ crossing time| [s]","Position RMS [m]"],
        [[label+' dt/2 + displacement/2',number(value['assignment_agreement']),number(value['f_A_difference']),number(value['A_over_initial_difference']),number(value['B_over_initial_difference']),number(value['mean_absolute_crossing_time_difference_s']),number(value['final_or_crossing_position_RMS_difference_m'])] for label,value in conv['timestep'].items()]+
        [["Classification tolerance",number(conv['classification']['assignment_agreement']),number(conv['classification']['f_A_difference']),number(conv['classification']['A_over_initial_difference']),number(conv['classification']['B_over_initial_difference']),number(conv['classification']['mean_absolute_crossing_time_difference_s']),number(conv['classification']['final_or_crossing_position_RMS_difference_m'])]])+
        "\n\nPosition difference compares final active or first-exit positions; it is not automatically a common-time displacement error. "
        f"Observed fraction uncertainties: `{uncertainties}`; minimum lower paired-effect CI: {number(lower)}. "
        "The displacement limit is also halved (5 R0 to 2.5 R0), including a B2 primary+secondary run. This tightens the control responsible for the effective microsteps; reducing only the nominal dt would miss that issue. "
        "One-seed tests quantify observed sensitivity, not a proven worst-case bound over all seeds. Per-ID branch-switch matrices, initial-population arrival differences and early clustering differences are retained in the numerical summary. "
        "The early-cluster relative L2 diagnostic uses a 1e-30 denominator guard: when the reference series is identically zero and the candidate is nonzero, its large stored value flags a failed relative gate and must not be interpreted as a meaningful percentage; use the retained absolute difference and both time series.\n\n"+
        table(["Timestep comparison","Original steps","Stricter steps","Original mean dt [µs]","Stricter mean dt [µs]"],
            [[label,value["effective_timestep"]["reference"]["accepted_steps"],value["effective_timestep"]["candidate"]["accepted_steps"],
              number(1e6*value["effective_timestep"]["reference"]["mean_effective_dt_s"]),
              number(1e6*value["effective_timestep"]["candidate"]["mean_effective_dt_s"])] for label,value in conv["timestep"].items()])+
        "\n\n"+
        table(["Alternate uniform-area case","A","B","U","f_A"],[[r['case'],r['N_A'],r['N_B'],r['N_unassigned'],number(r['f_A'])] for r in conv['sampling']])+
        "\n\n"+conv.get('population_size_scope',"Optional N512 results are recorded in the convergence summary."))
    grid = conv['acoustic_grid']
    grid_matrix = grid['branch_comparison']['branch_switch_matrix']
    grid_counts = [
        ["Coarse"]+[sum(grid_matrix[origin].values()) for origin in outcomes],
        ["Fine"]+[sum(grid_matrix[origin][destination] for origin in outcomes) for destination in outcomes],
    ]
    section("Acoustic-grid sensitivity",f"B2 fine direct k-Wave solve keeps the physical domain and controls, with spacing factor 0.8 and grid-count factor 1.25. "
        f"At {grid['common_sample_count']} common actual trajectory coordinates, relative L2 errors (fine norm denominator): pressure {number(grid['pressure_relative_L2'])}, "
        f"complex gradient {number(grid['gradient_relative_L2'])}, primary force {number(grid['primary_force_relative_L2'])}. "
        f"Fine-trajectory Δ f_A: {number(grid['branch_comparison']['f_A_difference'])}; assignment agreement {number(grid['branch_comparison']['assignment_agreement'])}. "
        f"Changes in initial-population arrival fractions: Δ A/N0 = {number(grid['branch_comparison']['A_over_initial_difference'])}, Δ B/N0 = {number(grid['branch_comparison']['B_over_initial_difference'])}. "
        "Force sensitivity can exceed pressure sensitivity. Two grids diagnose a tested numerical effect; they do not prove grid independence or certify a patient-scale acoustic model.\n\n"+
        table(["Acoustic grid","A","B","U","W"],grid_counts)+"\n\n"+
        table(["Coarse outcome / fine outcome","A","B","U","W"],
            [[label]+[grid_matrix[origin][destination] for destination in outcomes]
             for label,origin in zip(["A","B","U","W"],outcomes,strict=True)])+
        "\n\nA/B-to-U differences indicate an observed exit versus censoring at 6 s; they are not observed switches to the opposite branch. Equal conditional f_A can coexist with different absolute arrival counts. "+
        f"There are {grid['branch_comparison']['shared_outlet_crossings']} shared exits, with mean |Δ crossing time| {number(grid['branch_comparison']['mean_absolute_crossing_time_difference_s'])} s.\n\n"+
        link(str(Path(conv['run'])/'GRID_COMMON_POINTS.npz'),"Common-point complex fields and forces"))
    history_transitions = {origin:{destination:sum(row['production_branch']==origin and row['history_reference_branch']==destination
        for row in history['rows']) for destination in outcomes} for origin in outcomes}
    assert sum(sum(row.values()) for row in history_transitions.values()) == history['selected_count']
    section("History-force branch robustness",f"Selected {history['selected_count']} critical B2 particles automatically: prioritize B0/B2 assignment switches, then smallest |y| at first bifurcation crossing, then global ID. "
        "Only assigned B2 particles are eligible; this is a deliberately difficult, nonrepresentative sample. "
        "The validation-only reference updates local flow, pressure and gradient along its own moving path. It combines a variable-radius periodic APECSS translation closure in transformed time with an equilibrium-radius slow added-mass/Abel kernel. "
        "It is an **adiabatic two-scale approximation**, not fully resolved moving variable-radius history. Constant-radius tests, spectral refinement, independent closure holdouts and moving dt refinement constrain its numerical error.\n\n"
        "The conservative added momentum and no-slip history kernel follow [Garbin et al.](https://pages.jh.edu/aprospe1/publications/PapersPublished/GasBubbles/GarbinHistory.pdf) and "
        "[Cattaneo & Supponen, Eqs. 11 and 14](https://arxiv.org/html/2506.16413v2). Startup, fast/slow coupling, fluid-acceleration terms and wall hydrodynamics remain unresolved.\n\n"+
        table(["Comparison","Branch agreement"],[["Production vs reference C",f"{history['branch_agreement_rate']:.1%}"],
            ["C dt vs dt/2",f"{history['reference_dt_assignment_agreement']:.1%}"],
            ["Reduced-population production vs full population",f"{history['reduced_population_control_agreement']:.1%}"],
            ["Reference B (inertia) vs C (history)",f"{history['B_vs_C_history_only_branch_agreement']:.1%}"]])+
        "\n\n"+table(["Production outcome / reference C outcome","A","B","U","W"],
            [[label]+[history_transitions[origin][destination] for destination in outcomes]
             for label,origin in zip(["A","B","U","W"],outcomes,strict=True)])+
        "\n\nA/B denote observed first outlet crossings; U means no classified exit by the observation cutoff; W denotes wall loss. "
        "An A/B-to-U difference is an observed-versus-censored outcome difference, not an observed switch to the opposite branch. "
        "Agreement includes all selected particles and all four outcome labels; the transition counts distinguish these different causes of disagreement. "+
        f"\n\nReport classification: **{history['classification']}**. Production-vs-C includes inertia, microcycle averaging and different radial-table/periodic-closure interpolation errors. B-vs-C tests history within the same reference construction; reduced-population control identifies omitted geometric-contact confounding. "
        "Full per-ID crossing positions/times and displacement differences: "+link(str(Path(history['run'])/'SUMMARY.json'),"history summary")+". Production source has no history-module imports.")
    section("Performance",table(["Primary case","Median runtime [s]","Maximum runtime [s]"],[[label,number(row['median_elapsed_s'],6),number(row['maximum_elapsed_s'],6)] for label,row in perf.items()])+
        "\n\n"+link(str((out/'PERFORMANCE.csv').relative_to(ROOT)),"All per-case acoustic-loading/response/pair/flow/contact/integration/output timings")+
        ". Integration subtimers contain nested evaluations and must not be added as disjoint components. Concurrent cases share CPU quota. "
        "Profiling identified repeated mesh wall-distance evaluations during geometric contact projection. Exact per-row coordinate cache avoids duplicate VTK queries without quantization. Full B0 and short B1 trajectories before/after the cache were bitwise identical. "
        "A failed B1 contact step at budget 100 was replayed at unchanged dt and tolerance: budgets 300, 1,000 and 3,000 converged to identical results. Formal runs use budget 3,000; no forced acceptance or relaxed contact tolerance.\n\n"
        "A separate 0.04 s secondary-ON profile required 1,755 accepted substeps and 63,952 geometric pair-contact events. Repeated wall queries dominated; pair-force summation was only about 2.3% of the instrumented runtime. "
        "A diagnostic-only conservative interior-distance certificate reduced the measured runtime from 114.8 to 78.0 s while all 13 saved trajectory datasets, the contact-event stream and validity summary remained identical. "
        "The initial whole-batch certificate was not deployed. A later per-particle certificate retains the exact original wall correction for every uncertain point. "
        "B2/B3 replays from saved t=0.24 s states each preserved all 64 accepted steps, velocities, contacts and validity records bitwise, with measured runtimes reduced from 13.07/14.99 s to 6.07/6.45 s. "
        "This verified option is identified by a helper source hash in each enabled run. It skips provably inactive wall corrections without returning approximate distances or changing any force. "
        "The earlier pre-I.5 development batch `primary_secondary/20261003T115216806643Z` was stopped while incomplete to enable that exact optimization; its partial results are preserved and excluded. The replacement certified-wall batch `20261003T132343245938Z` restarted from the original initial conditions and became the unique full gold. Both gold workers completed naturally before I.5 began; this historical optimization is not a restart of the sealed baseline. "
        "The installed-version distance-query behavior is documented in the [VTK 9.7.1 implementation](https://raw.githubusercontent.com/Kitware/VTK/v9.7.1/Filters/Core/vtkImplicitPolyDataDistance.cxx). "
        "Evidence: "+link("results/stage4/performance/profile_B2_secondary/20261003T115538396917Z/SUMMARY.json","secondary profile")+"; "+link("results/stage4/performance/certified_wall_trial/20261003T120728567547Z/SUMMARY.json","early equivalence trial")+"; "+
        link("results/stage4/performance/later_particle_geometry_certificate/20261003T131349331620Z/SUMMARY.json","later per-particle equivalence trial")+".\n\n"+
        "A separate component-wise contact NNLS optimization was rejected. Across all 3,519 shadow-checked queries, velocity differences were at most 8.95e−16 m/s (relative L2 at most 2.76e−15), yet the maximum coordinate difference at the same saved times reached 49.1 µm by 0.04 s. "
        "Contact/validity histories also differed. The original global contact solve remains in production. This establishes individual-trajectory sensitivity in the tested secondary-ON diagnostic; final branch-fraction sensitivity was not measured by that short test. "+
        link("results/stage4/performance/component_contact_shadow/20261003T130544461050Z/SUMMARY.json","contact roundoff diagnostic"))
    profile=stage("performance_refactor/profiles_report")
    limiter=stage("performance_refactor/limiter_analysis")
    gold_folder=Path("results/stage4/performance_refactor/gold_baseline")
    gold_rows=[]
    for label in ["B2","B3"]:
        row=read(gold_folder/"cases"/label/"particles/SUMMARY.json")
        b=row["branch_metrics"]
        gold_rows.append([label,row["accepted_steps"],number(row["duration_s"]/row["accepted_steps"]*1e6),
            number(row["elapsed_s"]/3600),b["N_A"],b["N_B"],b["N_unassigned"]])
    section("Performance bottleneck discovered",table(["Gold case","Accepted steps","Mean dt [µs]","Wall time [h]","A","B","Unassigned"],gold_rows)+
        "\n\nThe 0.5 ms nominal clock has 12,000 macro steps over 6 s, but the global displacement minimum produces roughly 0.9 million accepted steps. "
        "The minimum dt can also be a macro-boundary remainder; the limiter audit distinguishes these cases. Contact position projection, global NNLS and wall work dominate different segments; pair-force acceleration is not the main opportunity. "
        "The updated P0 profile includes the exact interior certificate already used in gold, so earlier uncertified timings are not used as the I.5 speedup denominator.\n\n"+
        table(["Segment","Accepted steps","Limiter IDs covering 90%","Median stiff / active","Median contact component"],
            [[r["segment"]["id"],r["accepted_steps"],r["particles_covering_90_percent"],
              f"{r['median_particles_requiring_substeps']:.0f} / {r['median_active_particles']:.0f}",r["median_limiter_contact_component"]] for r in limiter["segments"]])+
        "\n\n"+link(str(Path(profile["run"])/"PROFILE_STAGE4I5_BASELINE.md"),"Full exclusive timer and cProfile report")+"; "+
        link(str(Path(limiter["run"])/"TIMESTEP_LIMITER_ANALYSIS.md"),"Every-step limiter analysis")+".\n\n"+
        "Near/contact pair fractions weighted by evaluation count, accepted physical time and wall time are different quantities. A step can simultaneously contain far-field and near-contact pairs. The fraction of runtime with any near/contact pair is reported separately from weighted pair fractions. "+
        link(str(Path(limiter["run"])/"PERFORMANCE_VS_PHYSICS_VALIDITY.md"),"Performance versus model-validity evidence")+
        ". Advancing an approximate far-field force into close configurations creates numerical stiffness; more computation and smaller dt do not establish higher physical accuracy.\n\n"+
        table(["Segment","Instrumented step time with any near/contact pair [%]","Time-weighted near pair fraction [%]","Time-weighted contact pair fraction [%]"],
            [[r["segment"]["id"],number(100*r["validity_weighting"]["measured_step_wall_time"]["any_NEAR_CONTACT_or_CONTACT"]),
              number(100*r["validity_weighting"]["accepted_physical_time"]["NEAR_CONTACT"]),
              number(100*r["validity_weighting"]["accepted_physical_time"]["CONTACT"])] for r in limiter["segments"]])+
        "\n\nThe first percentage measures runtime coverage by at least one near/contact pair; the last two average population pair fractions over accepted physical time. These columns do not partition runtime and are not a full-six-second instrumented census.")
    section("Stage 4I.5 acceleration",f"Recommended J/K backend: **{decision['recommended_backend_for_JK']}**. "
        "P0 retains the original global RK2/NNLS methods. P1 tests exact-coordinate wall caching, a persistent vtkStaticCellLocator and a conservative cKDTree skin list separately and in combination. "
        "The static locator is built once and reuses query objects; the old VTK implementation already built its locator once, so no nonexistent rebuild bottleneck is claimed. "
        "Neighbor broad-phase reuse keeps the original exact narrow-phase threshold and sorted order; stress tests require zero missed contacts. "
        "The additional ordered CPU contact sweep compiles only the original sequential correction loop; it reuses the original NumPy norm operation and keeps the global NNLS, candidate ordering and physical tolerances. Compiler, binary and BLAS identities are recorded. "
        "Warp is isolated and optional, with closed-mesh/winding, 10,000-point, transfer-cost and trajectory checks. The default environment remains usable without it.\n\n"+
        table(["Candidate","Completed matched segments","Strict correctness","Component-timed aggregate speedup"],
            [[r["backend"],r["complete_segments"],r["strict_correctness"],number(r["aggregate_speedup"])] for r in decision["candidates"]])+
        "\n\nThe component-timed runs include millions of per-point VTK wrappers for the original backend; instrumentation cost can depend on the query API. Promising strict candidates therefore also undergo the same four short windows with component wrappers removed before adoption. Common per-step recording, timestep and classifier timers remain, and timer-removal trajectories are checked against their instrumented counterparts. These timings are not full-production runtime measurements.\n\n"+
        table(["Candidate","Strict minimal-timer confirmation","Confirmed aggregate speedup"],
            [[r["backend"],r["strict_correctness"],number(r["aggregate_speedup"])]
             for r in (decision.get("minimal_timer_confirmation") or {}).get("candidates",[])])+
        "\n\nIndependent 10,000-point wall fixture (query speed includes host adapter transfers; decisions retain the original tolerances):\n\n"+
        table(["Wall backend","Geometry gate","Query speedup","Max distance error [m]","Max normal angle [deg]","Projection decision mismatches"],
            [[name,r["pass_geometry"],number(r["query_speedup"]),number(r["max_distance_error_m"]),
              number(r["max_normal_angle_deg"]),r["projection_decision_mismatches"]]
             for name,r in decision["independent_prerequisite_gates"].items() if name in ["vtk","warp"]])+
        "\n\n"+decision["p3"]["scope"]+" An SDF is not substituted for exact geometry. The 2× target is an engineering aim, not a reason to relax correctness. "
        "No force/velocity cap, cutoff, pressure change, radius change or relaxed contact/wall tolerance was introduced.\n\n"+
        link(str(Path(decision["run"])/"STAGE4I5_BACKEND_DECISION.md"),"Backend decision and rejected candidates")+"; "+
        link("docs/ASSIGNED_PARTICLE_DEACTIVATION_AUDIT.md","Assigned-particle boundary audit")+
        ". Existing absorbing classification boundaries are preserved for comparison, but physical terminality is unproven. No additional deactivation or plateau-based stop is enabled. J validates its declared finite-observation objective; unassigned particles prevent interpreting a plateau as endpoint completion.")
    regressions=[r for candidate in decision["candidates"] for r in candidate["regression"]]
    section("Gold-baseline regression",table(["Backend / segment","Max position error [m]","Branch agreement","New crossings (reference)","Max new crossing-time error [s]","Exact contact sequence","Speedup"],
        [[r["backend"]+" / "+r["segment"],number(r["max_position_error_m"]),number(r["branch_agreement_all_ids"]),
          r["new_crossings_reference"],number(r.get("maximum_new_crossing_time_error_s")),r["contact_sequence_identical"],number(r["speedup"])] for r in regressions])+
        "\n\nBoth original Stage 4I cases naturally completed before any refactor. Their source/config/input hashes, trajectories, crossings and event logs are sealed read-only and independently backed up. "
        "The original long baseline was never rerun. P0 short replay endpoints are checked directly against stored gold, and candidate state/contact comparisons use those matched replays. "
        "A segment with no new outlet event cannot establish long-horizon branch robustness merely by reporting 100% assignment agreement. "
        "No full optimized P3 case was run when local-subcycling eligibility was absent; no unmeasured full-run optimized agreement or runtime is invented. "
        "Original GPU time, complete rejected-step/max-dt instrumentation and some geometric counters were not recorded and remain explicitly missing in the manifest.\n\n"+
        link(str(gold_folder/"GOLD_BASELINE_MANIFEST.json"),"Immutable gold manifest")+"; "+
        link(str(gold_folder/"STAGE4I_BASELINE_SUMMARY.md"),"Gold metrics and instrumentation limits"))
    validity = []
    for row in stages['primary_only']['rows']+stages['primary_secondary']['rows']:
        v = read(Path(row['run'])/'physics_validity.json')
        s = read(Path(row['run'])/'particles/SUMMARY.json')
        validity.append([row['case']+(' + secondary' if row['secondary'] else ' primary'),v['validity_status'],v['physics_status'],s['wall_hit_particle_count'],number(v['invalid_accepted_duration_s']),number(v['pair_regime']['near_field_fraction']),number(v['pair_regime']['contact_fraction'])])
    section("Model validity",table(["Run","Validity banner","Physics scope","Particles with wall projection","Any-reason flagged time [s]","Near-field pair fraction","Contact pair fraction"],validity)+
        "\n\nFlagged time means an accepted endpoint had at least one configured validity reason, including any near-wall particle; it is not the fraction of invalid individual trajectories. "
        "The substantial wall-projection population directly limits physical branch-distribution claims. Neither history agreement under geometric walls nor a smooth trajectory resolves missing wall hydrodynamics. \n\n"+
        table(["Run","Maximum saved constrained-slip Re","Saved active samples with Re ≥ 1"],
              [[r['case']+(' + secondary' if r['secondary'] else ' primary'),number(r['saved_slip_reynolds']['measurements']['constrained']['maximum']),
                r['saved_slip_reynolds']['measurements']['constrained'].get('sample_count_Re_at_least_1',0)] for r in model_evidence])+
        "\n\nRe = 2ρR0|v − u_flow|/µ uses the velocity after geometric constraints. These are saved-state peaks, so unrecorded substep peaks can be larger. Sample counts are neither particle counts nor time fractions. Free velocities without contact reactions are recorded separately in the machine-readable diagnostics. No finite-Re drag correction has been added. "+
        "\n\nProduction remains linear one-way acoustics, adaptive radial response, cycle-averaged primary Bjerknes, classical secondary Bjerknes, quasi-steady Stokes transport, imported external flow and geometric nonpenetration. "
        "Pair fractions are evaluation-pair weighted including RK stages; geometric separation labels are not a complete validation of large-separation retardation or near-wall physics. "
        "No streaming, history in production, fusion, adhesion, lubrication, nonspherical deformation, electromechanical FEM, skull/brain acoustics, acoustic feedback or bubble-flow feedback was added. The concentration surrogate is normalized finite-ensemble occupancy in the nominal 1–3 mm daughter ROI, not a clinical particles/mL prediction.")
    section("What is now safe to claim",("The tested phase/amplitude controls change the simulator's conditional branch preference with paired-seed evidence. " if stats['statistical_directional_control_supported'] else "The current tested controls do not establish the requested directional conditional branch enrichment with the required paired-seed evidence. ")+
        "The existing particle physics uses the steady SimVascular velocity field through the unchanged FlowField API. "
        "Any claim of numerical/model robustness is restricted to the component statuses above. Conditional exit fractions, target arrival counts and unassigned fraction must be reported together. "
        "No claim is made about the real percentage of SonoVue reaching a target cerebral vessel.")
    section("What remains uncertain","CFD mesh/local near-wall accuracy; acoustic-grid convergence beyond two resolutions; wall hydrodynamics; near-contact/classical pair physics; slow/fast history coupling and unmodelled fluid acceleration; selected-reference generalization; longer observation and release of trapped bubbles; particle-size/frequency distributions outside the narrow extended table; isolated phase-only causality; injection concentration; experimental calibration. "
        "B1 symmetry and B2/B3 mirror evidence are evaluated through the full paired statistics rather than a single seed. Zero mass imbalance, smooth visualization and GPU consistency are not proofs of local physical accuracy.")
    section("Readiness for real SimVascular vascular geometry / MeVO",("The simplified coupling bridge passes the stated scoped numerical/model gates. " if success else "The simplified coupling bridge executes, but at least one control/robustness gate remains a warning; quantitative branch-control validation is not established. ")+
        "The input/export adapters and run identities can support later geometry work. No patient data, asymmetric vascular network or pulsatile simulation was started in this stage.")
    section("Recommended Stage 5","Before a new geometry study, resolve the warning components in STAGE4_ACCEPTANCE, preserve observation censoring, and choose control amplitudes/windows with interpretable arrival statistics. "
        "A separately authorized next stage may investigate an asymmetric nonpatient SimVascular network, mesh/acoustic convergence and then pulsatile flow. Experimental and patient-specific claims require additional calibration. These are recommendations only; Stage 5 has not been executed.\n\n"
        "The eight scientific figure captions and exact sources are in "+link(str((out/'figures/FIGURE_MANIFEST.json').relative_to(ROOT)),"figure manifest")+"; six additional performance figures are in "+link(str((out/'performance_figures/FIGURE_MANIFEST.json').relative_to(ROOT)),"performance figure manifest")+". ParaView active-particle series, native flow, labelled vessel and array exports are in "+link(str((out/'paraview/EXPORT_MANIFEST.json').relative_to(ROOT)),"ParaView manifest")+".")
    report = "# Stage 4 — Y-Bifurcation PZT Branch Steering & SimVascular Coupling\n\n"+"\n".join(sections)
    (out / "STAGE4_VALIDATION_REPORT.md").write_text(report.replace("](../", "](../../../../"))
    (ROOT / "docs/STAGE4_VALIDATION_REPORT.md").write_text(report)
    cards = "\n".join(f'<article><img loading="lazy" src="{Path(f["png"]).relative_to(out.relative_to(ROOT))}" alt="{escape(f["name"])}"><p>{escape(f["caption"])}</p></article>' for f in figures)
    viewer = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Stage 4 | Y-bifurcation control evidence</title><style>
body{background:#05090f;color:#ecf3ff;font:17px/1.65 system-ui,sans-serif;max-width:1450px;margin:auto;padding:30px}
a{color:#59dbca}h1{line-height:1.2}header{max-width:1100px}article{margin:50px 0}img{width:100%;height:auto}p{color:#b9cddd}
code{color:#eeb770;overflow-wrap:anywhere}nav{display:flex;gap:24px;flex-wrap:wrap}</style><header>
<h1>Stage 4 — Y-bifurcation PZT branch steering</h1><p>256 bubbles × 10 paired seeds · 6 s observation · steady SimVascular flow</p>
<p>Scientific state: <code>"""+escape(acceptance["stage_state"])+"""</code></p>
<p>Conditional exit preference, target arrival count and unassigned fraction are separate outcomes. Controls vary phase and amplitude.
Geometric contact is not fusion; the history reference is validation-only and uses an adiabatic two-scale approximation.</p>
<nav><a href="STAGE4_VALIDATION_REPORT.md">Validation report</a><a href="STAGE4_ACCEPTANCE.json">Acceptance evidence</a>
<a href="SECONDARY_COMPARISON.json">Secondary metrics</a><a href="PERFORMANCE.csv">Runtime table</a>
<a href="paraview/EXPORT_MANIFEST.json">ParaView source manifest</a></nav></header>"""+cards+"</html>\n"
    (out/"OPEN_RESULTS.html").write_text(viewer)
    readme = ROOT/"README.md"
    body = readme.read_text()
    begin = body.find("## Stage 4 — Y-bifurcation branch steering")
    end = body.find("## Stage 3 — minimal validated simulator",begin)
    readme_updated = begin>=0 and end>begin
    if readme_updated:
        section_text = ("## Stage 4 — Y-bifurcation branch steering\n\n"
            f"状态：`{acceptance['stage_state']}`（`{acceptance['status']}`）。计算、适用性警告及各项证据见最终报告。"
            "条件分支偏向、实际目标到达量和观察截止时未出流比例须同时报告；近壁水动力学与临床输送效果尚未验证。\n\n"
            "- [Stage 4 验证报告](docs/STAGE4_VALIDATION_REPORT.md)与 [验收证据](STAGE4_ACCEPTANCE.json)\n"
            f"- [科学结果图、六张性能图与 ParaView 数据入口]({(out/'OPEN_RESULTS.html').relative_to(ROOT)})\n"
            "- [运行与接续说明](docs/STAGE4_QUICKSTART.md)、[预先声明的执行范围](docs/STAGE4_PLAN.md)\n"
            "- [2026-10-03 历史进度与运行诊断](docs/STAGE4_RUNTIME_NOTE_20261003.md)\n\n")
        readme.write_text(body[:begin]+section_text+body[end:])
    # Preserve exact evidence pointers independently of the mutable latest-state files.
    write(out / "STAGE_RUNS.json",stages)
    finish(out,dict(status=acceptance['status'],stage_state=acceptance['stage_state'],figures=figures,
                    paraview_collections=collections,report="docs/STAGE4_VALIDATION_REPORT.md",acceptance="STAGE4_ACCEPTANCE.json",
                    viewer=str((out/"OPEN_RESULTS.html").relative_to(ROOT)),
                    readme_updated=readme_updated,
                    report_sha256=sha(ROOT/"docs/STAGE4_VALIDATION_REPORT.md"),acceptance_sha256=sha(ROOT/"STAGE4_ACCEPTANCE.json")))


if __name__=="__main__":
    main()
