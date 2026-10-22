"""Deterministic sparse 6N affine resistance assembly, stationary WALL, Fnh=0."""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from .hydrodynamic_resistance import sphere_self_diagonal, evaluate_block


@dataclass
class ResistanceSystem:
    ids: tuple
    matrix: sparse.csr_matrix
    rhs: np.ndarray
    self_diagonal: np.ndarray
    free: np.ndarray
    blocks: list
    block_rows: list

    def dissipation(self, velocity):
        u = np.asarray(velocity)
        slip = u-self.free
        result = dict(self=float(np.dot(slip*self.self_diagonal, slip)), wall=0., pair=0.)
        for record, row in zip(self.blocks, self.block_rows):
            result[record['kind'].lower()] += float(record['coefficient_kg_s']*(row@u)**2)
        result['total'] = sum(result.values())
        return result


def assemble_resistance_system(particles, free_velocities, wall_blocks=(), pair_blocks=(), *, mu):
    ids = tuple(sorted(particles))
    if not ids:
        raise ValueError('At least one particle required')
    index = {i: 6*k for k, i in enumerate(ids)}
    diag = np.concatenate([sphere_self_diagonal(particles[i], mu) for i in ids])
    free = np.concatenate([np.asarray(free_velocities[i], dtype=float) for i in ids])
    if free.shape != (6*len(ids),) or not np.isfinite(free).all():
        raise ValueError('Each free generalized velocity must have six finite components')
    matrix = sparse.diags(diag, format='csr')
    blocks, rows = [], []
    specs = list(wall_blocks)+list(pair_blocks)
    # Reject foreign objects before sorting/accessing attributes.
    from .hydrodynamic_resistance import PhysicalNearField
    if any(type(s) is not PhysicalNearField for s in specs):
        raise TypeError('VALIDATION_BLOCK_CANNOT_ENTER_PHYSICAL_SOLVER')
    specs.sort(key=lambda s: (s.particle_j_id is not None, min(s.particle_i_id, s.particle_j_id if s.particle_j_id is not None else s.particle_i_id),
                              max(s.particle_i_id, s.particle_j_id if s.particle_j_id is not None else s.particle_i_id), tuple(s.normal)))
    for spec in specs:
        zeta, normal, record = evaluate_block(spec, particles, mu)
        row = np.zeros(len(diag))
        i = index[spec.particle_i_id]; row[i:i+3] = normal
        if spec.particle_j_id is not None:
            j = index[spec.particle_j_id]; row[j:j+3] = -normal
        if zeta:
            r = sparse.csr_matrix(row.reshape(1, -1))
            matrix += zeta*(r.T@r)
        record['normal'] = normal.tolist()
        blocks.append(record); rows.append(row)
    return ResistanceSystem(ids, matrix, diag*free, diag, free, blocks, rows)
