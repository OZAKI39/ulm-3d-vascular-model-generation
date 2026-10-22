"""Explicit source evidence gate and pressure transfer between physical planes."""
import numpy as np


class SourceAmbiguous(ValueError):
    pass


def require_source(source_evidence):
    """A serialization root alone is not an identified hydraulic source."""
    if (source_evidence.get('hydraulic_source_id') is None
            or source_evidence.get('identification_kind') not in {'DATASET_BOUNDARY_ANNOTATION', 'EXPLICIT_MODEL_ASSUMPTION'}
            or not source_evidence.get('evidence_reference')):
        raise SourceAmbiguous('A_NETWORK_SOURCE_AMBIGUOUS')
    return int(source_evidence['hydraulic_source_id'])


def cap_pressures(real_pressure, outward_flow, extension_resistance):
    """For outward flow, the distal artificial cap is at p_real - R_ext Q_out."""
    p, q, r = np.broadcast_arrays(real_pressure, outward_flow, extension_resistance)
    if not all(np.all(np.isfinite(x)) for x in (p, q, r)) or np.any(r < 0):
        raise ValueError('Invalid pressure transfer inputs')
    return p-r*q


def shift_pressure_gauge(pressure):
    p = np.asarray(pressure, dtype=float)
    if p.size == 0 or not np.all(np.isfinite(p)):
        raise ValueError('Finite pressure set required')
    shifted = p-p.min()
    if not np.allclose(p[:, None]-p, shifted[:, None]-shifted, rtol=1e-12, atol=1e-12):
        raise ValueError('Gauge shift changed pressure differences')
    return shifted
