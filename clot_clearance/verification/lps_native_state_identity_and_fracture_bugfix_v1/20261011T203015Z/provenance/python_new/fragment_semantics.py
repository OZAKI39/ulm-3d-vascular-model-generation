"""Versioned attachment-aware component metrics; legacy fields remain unchanged.

Denominator is total input material mass. Callers must exclude wall anchors.
Tiny (1–2 node) components are unresolved numerical debris, not microemboli.
"""
import numpy as np

SEMANTICS_VERSION = 'ATTACHMENT_AWARE_V2'
LEGACY_FIELDS = {
    'debris_fraction': 'all component mass fraction with 1–2 nodes (equal-volume legacy)',
    'debris_1_2node_fraction': 'all component mass fraction with 1–2 nodes',
    'resolved_fragment_count': 'all components with at least 3 nodes, including attached',
    'largest_free_fraction': 'LEGACY_MISNOMER: largest ALL component / total material',
    'FI_free': 'LEGACY_MISNOMER: 1 - largest ALL component / total material',
}

def component_metrics(counts, attached, dx, rho=1000., *, component_mass=None, component_volume=None):
    counts=np.asarray(counts,dtype=int);attached=np.asarray(attached,dtype=bool)
    if counts.ndim!=1 or counts.shape!=attached.shape or np.any(counts<=0):
        raise ValueError('Positive component node counts and matching attachment flags required')
    volume=counts*dx**3 if component_volume is None else np.asarray(component_volume,float)
    mass=rho*volume if component_mass is None else np.asarray(component_mass,float)
    if volume.shape!=counts.shape or mass.shape!=counts.shape or np.any(volume<=0) or np.any(mass<=0):
        raise ValueError('Positive component masses/volumes must match components')
    total=mass.sum();free=~attached;tiny=counts<=2
    frac=lambda mask:float(mass[mask].sum()/total)
    largest=lambda mask:float(mass[mask].max(initial=0)/total)
    row=dict(fragment_semantics_version=SEMANTICS_VERSION,
        all_tiny_component_mass_fraction=frac(tiny),
        detached_tiny_component_mass_fraction=frac(tiny&free),
        attached_tiny_component_mass_fraction=frac(tiny&attached),
        detached_component_count=int(free.sum()),
        resolved_detached_component_count=int((free&~tiny).sum()),
        largest_detached_nodes=int(counts[free].max(initial=0)),
        largest_detached_mass_fraction=largest(free),
        largest_attached_component_mass_fraction=largest(attached),
        detached_mass_fraction=frac(free),detached_mass_kg=float(mass[free].sum()),
        attached_mass_kg=float(mass[attached].sum()),total_material_mass_kg=float(total),
        material_node_count=int(counts.sum()),
        median_detached_component_nodes=float(np.median(counts[free])) if free.any() else None,
        legacy_fields_semantics='ALL_COMPONENTS_INCLUDING_ATTACHED_V1')
    for diameter in [200,300]:
        # Volume comparison avoids an equality-rounding error in cube roots.
        selected=free&(volume >= np.pi*(diameter*1e-6)**3/6)
        row['N_detached_ge'+str(diameter)]=int(selected.sum())
        row['mass_fraction_detached_ge'+str(diameter)]=frac(selected)
    return row

def classify_detached_v2(m):
    """New opt-in labels. Never reinterpret/overwrite legacy failure_mode or AI labels."""
    if m['detached_tiny_component_mass_fraction']>=.5:return 'DETACHED_UNRESOLVED_DEBRIS_DOMINANT'
    if m['largest_detached_mass_fraction']>=.8:return 'WHOLE_CLOT_DETACHMENT'
    if m['resolved_detached_component_count']>0:return 'DETACHED_RESOLVED_COMPONENTS'
    if m['detached_component_count']>0:return 'ATTACHED_CORE_WITH_LOCAL_DEBRIS_SHEDDING'
    if m.get('bulk_broken_fraction',0)>0 or m.get('interface_broken_fraction',0)>0:return 'ATTACHED_DAMAGE_WITHOUT_DETACHMENT'
    return 'NO_FAILURE'
