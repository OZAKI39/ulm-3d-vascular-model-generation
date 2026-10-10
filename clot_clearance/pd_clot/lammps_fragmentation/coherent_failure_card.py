"""Liu 2024 only: distinguish measured quantities from inverse-model ranges."""
from pathlib import Path
import math
from .config import BASELINE, ROOT

OUTPUT = ROOT / 'verification/lammps_coherent_failure_force_time_map_v1'
RESULTS = ROOT / 'results/lammps_coherent_failure_force_time_map_v1'
DOI = '10.1016/j.jmps.2024.105858'

def failure_card():
    return dict(card_id='COHERENT_FAILURE_CARD_LIU2024_V1',source_DOI=DOI,
        sample='Recalcified bovine whole blood; normal RBC/platelet condition',
        preparation='CPDA-1 blood; blood:recalcification solution 10:1; final 22.5 mM NaCl and 30 mM CaCl2; 37 C 60 min; 4 C overnight; tests next day',
        source_rate_s_inv=.01,bulk_loading_speed_mm_s=.04,peeling_loading_speed_mm_s=.02,
        bulk=dict(measured_energy_J_m2=18.27,measured_energy_SD_J_m2=1.45,observed_strength_Pa=4000.,
                  cohesive_strength_range_Pa=[4000.,8000.],intrinsic_work_range_J_m2=[6.,8.],
                  identity='Energy and observed strength experimental; strength interval and intrinsic work constrained FEM inference',
                  intrinsic_work_condition='Paper also invokes S_bulk >5.2 kPa when narrowing intrinsic work to 6-8 J/m2'),
        interface=dict(production_reference_substrate='PORCINE_SKIN_PEELING_MODEL_PROXY',observed_peeling_strength_Pa=1700.,
                       cohesive_strength_range_Pa=[1700.,3500.],intrinsic_work_range_J_m2=[1.3,1.6],
                       identity='Strength/work intervals inferred with porcine skin, Fig 8e; not calibrated aortic adhesion',
                       inner_aorta_separate_measurement=dict(mean_J_m2=1.35,SD_J_m2=.39,used_to_calibrate=False)),
        interface_failure_resistance_less_than_bulk=True,
        energies_used_as_critical_stretch=False,calibration_target='Macro nominal strength in declared small cells; no energy or curve fit',
        scope='EXPLORATORY_ONLY',biological_quantitative_validation='NOT_ESTABLISHED')

def elastic_mapping():
    return dict(status='ESTABLISHED',identity='LPS_ELASTIC_MAPPING_FROM_LIU2024',supplement='FOUND',
        actual_parameter_table='S2 (S1 is nomenclature)',C1_Pa=470.,C2_Pa=170.,C3_Pa=-12.,D_Pa_inv=.00006,
        G0_Pa=940.,K0_Pa=2/.00006,
        formula_source='Main Eq (3), I1 explicitly defined as first deviatoric invariant; undamaged eta=1',
        derivation='At J=1 simple shear gamma, I1bar-3=gamma^2, W=C1*gamma^2+O(gamma^4)=G0*gamma^2/2. At pure volume change, W=(J-1)^2/D=K0*(J-1)^2/2. Thus G0=2*C1 and K0=2/D.',
        units='1 kPa=1000 Pa; 0.06 kPa^-1=6e-5 Pa^-1',
        scope='Undamaged infinitesimal elastic mapping only; Yeoh nonlinearity and Mullins dissipation are not reproduced')

def configuration():
    c={**BASELINE,'K_Pa':2/.00006,'G_Pa':940.,'material_id':'LPS_ELASTIC_MAPPING_FROM_LIU2024',
       'load_provider':'CURRENT_POSITION_STREAMING_PROVIDER_V2','wall_interface_model':'BREAKABLE_INTERFACE_LIU2024_V1',
       'nominal_average_traction_Pa':None,'load_budget_N':None,'s00':None,
       'max_concurrent_cases':2,'mpi_ranks':2,'numerical_drag_rate_s_inv':500.,
       'pulse_rise_fraction':.25,'pulse_hold_fraction':.5,'pulse_decay_fraction':.25,
       'post_pulse_crossings':2.,'initial_amplitudes_Pa':[10.,100.,1000.],
       'duration_factors':[.5,2.,8.],'max_refinement_cases':12,'max_formal_cases':39,
       'source':'Liu 2024 strength-mapped LPS; analytic engineering load field',
       'viscoelastic_rate_dependence':'NOT_MODELED_V1'}
    c['T_cross_s']=c['radius_m']/math.sqrt(c['G_Pa']/c['rho_kg_m3'])
    return c
