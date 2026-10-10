"""Literature-informed elastic scales, explicitly assumed failure and density."""
import math

def prepare_material(request):
    r=request; V=r['dx_m']**3
    return dict(material_id='LPS_CLOT_LIKE_GULTEKIN_SCALE_V1',
        K=r['K_Pa'], G=r['G_Pa'], rho=r['rho_kg_m3'], s00=r['s00'], alpha=r['alpha'],
        delta=r['delta_m'], dx=r['dx_m'], nodal_volume_m3=V, nodal_mass_kg=r['rho_kg_m3']*V,
        contact_micromodulus_N_m6=18*r['K_Pa']/(math.pi*r['delta_m']**4),
        material_scope='EXPLORATORY_ONLY; linear state-based PD approximation, not equivalent to Ogden',
        K_source='Round A material card; assumed volumetric modulus in Gultekin et al. 2024',
        G_source='Round A Ogden small-strain shear scale mu=1040 Pa',
        density_source='ENGINEERING_DEMO_ASSUMPTION',
        critical_stretch_source='DEMO_CRITICAL_STRETCH_V1; not inferred from Gc',
        alpha_source='Engineering choice: constant critical stretch',
        source_DOI='10.1007/s13239-024-00733-2', patient_specific=False,
        biological_quantitative_validation='NOT_ESTABLISHED',
        SI_force_adapter='Multiply official internal force density by nodal volume before physical loading')
