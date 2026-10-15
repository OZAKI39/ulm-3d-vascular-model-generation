from common import *
import time,shutil
assert json.loads((O/'BASELINE_LPS_ENERGY_FORCE_CONSISTENCY.json').read_text())['status']=='PASS'
assert json.loads((O/'FIXED_DAMAGE_FORCE_GRADIENT_AUDIT.json').read_text())['status']=='PASS'
assert not (O/'DESIGN_FREEZE.json').exists()
freeze=dict(task=C['task'],run_id=C['run'],epoch=time.time(),candidate='ENDPOINT_STAR_VARIATIONAL_SERIES_COMPLIANCE_LPS_V1',candidate_count=1,
 K_Pa=K,G_Pa=G,rho_kg_m3=1000,s_onset=S0,alpha=0,dx_m=DX,horizon_m=DELTA,reference_sha256=sha(O/'provenance/reference.npz'),
 controls=[dict(name='instant_reference',s_complete=None),dict(name='soft_sc018',s_complete=.18),dict(name='soft_sc024',s_complete=.24)],
 path=dict(type='all-node controlled isotropic dilation x=(1+s)X',knots=[0,.15,.165,0,-.02,.165,.245],intervals=[60,30,40,10,50,160],event_points=[.150000001,.18,.24],purpose='kinematically compatible constitutive diagnostic, NOT equilibrium specimen or natural crack',length_scale_m=DX),
 auxiliary=dict(type='same 341-node two-half translation along X',plane_x_m=-30e-6,knots_over_dx=[0,.18,0,.18,.245],intervals=[30,20,20,25],s_complete=.24,area_m2=(7*DX)**2,purpose='nonuniform bond damage and observable partial planar separation, not natural fracture; no new specimen'),
 energy_update='At fixed old damage integrate -sum(F dot dx) over all controlled DOFs with 2-point Gauss; then update damage at fixed new coordinates; signed D=U_before-U_after; no abs or clipping. Splitting error assessed against exact uniform-dilation energy.',
 integration_order=2,physical_time=None,inertia=False,damping=False,contact=False,
 gates=dict(relative_work_balance=2e-8,uniform_exact_energy_relative=.012,negative_damage_J_tolerance=1e-20,damage_healing_tolerance=1e-14,zero_failed_edge_force_abs_N=1e-20),
 strain_scope=dict(calibrated_linear_strain_ceiling=None,diagnostic_max_bond_stretch=.245,statement='s0=.15 and 24.5% excursions are mathematical diagnostics, NOT evidence of small-strain biological validity'),
 high_energy_budget=dict(protocol='monotone two-half opening Delta in [0,.245 dx]',target_lower_J_m2=6.,method='Cauchy-Schwarz / Loewner upper bound on controlled work; no Gc mapping or sc extrapolation'),
 budget=dict(cpu_threads=1,max_wall_s=180,formal_path_cases=4,no_retries_for_unfavorable_science=True,no_native_fracture=True,no_V3=True,no_16000=True,no_GPU=True),
 sources={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')})
dump(O/'DESIGN_FREEZE.json',freeze)
dump(O/'RECOVERY_RECORD.json',dict(scientific_repeats=1,logging_failures=1,reason='Post-validation final console print used numpy.bool in standard json.dumps. All native/gradient audits and assertions had already completed and were saved. Corrected print conversion only; no numerical rerun.',original_log_sha256=sha(O/'checks.log'),recovery='Reused final-double-fix baseline twice: first draft exposed cross-component energetic coupling with shared internal elastic states; corrected allocation to independent endpoint stars within the same variational compliance mechanism, then reran necessary gates before any formal path. Draft and failure evidence retained.'))
print('FROZEN',sha(O/'DESIGN_FREEZE.json'))
