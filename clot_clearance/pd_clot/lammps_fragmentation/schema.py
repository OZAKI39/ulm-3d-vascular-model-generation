"""Input checks; future providers supply the same physical force table."""
import numpy as np

def validate(request):
    if request.get('scope') != 'EXPLORATORY_ONLY': raise ValueError('EXPLORATORY_ONLY required')
    if request.get('geometry_type') not in ('HEMISPHERE','STL','VTK','FEM_MESH'):
        raise ValueError('Unsupported geometry_type')
    for k in ('dx_m','delta_m','K_Pa','G_Pa','rho_kg_m3','s00','dt_s','ramp_time_s','load_budget_N'):
        if not np.isfinite(request[k]) or request[k] <= 0: raise ValueError(f'Invalid {k}')
    if request['delta_m'] <= request['dx_m']: raise ValueError('Horizon must exceed spacing')
    if request.get('alpha') != 0: raise ValueError('V1 freezes constant critical stretch: alpha=0')
    if request['hold_time_s'] < 0: raise ValueError('Negative hold time')
    if request.get('units') != dict(length='m',force='N',traction='Pa',mass='kg',time='s'):
        raise ValueError('Explicit SI units required')
    return dict(request)

def validate_force_table(ids, forces, node_count):
    ids=np.asarray(ids); forces=np.asarray(forces)
    if forces.shape != (len(ids),3) or not np.isfinite(forces).all(): raise ValueError('Invalid force array')
    if len(np.unique(ids)) != len(ids) or np.any(ids<1) or np.any(ids>node_count): raise ValueError('Invalid surface_node_id')
    return True
