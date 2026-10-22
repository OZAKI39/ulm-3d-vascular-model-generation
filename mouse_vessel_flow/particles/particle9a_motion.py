"""P6.5 accepted-time integration plus one sphere/nearest-wall affine block."""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from .particle_shapes import Sphere
from .resistance_assembly import ResistanceSystem
from .particle81_simulation import SavedTrajectoryStepper
from .particle81_cache import cached_step_to
from .planar_wall_hydrodynamics import planar_wall_affine_block
from .open_boundary_rim import rim_witness


@dataclass
class PlanarResistanceSystem(ResistanceSystem):
    planar_matrix: object=None
    planar_rhs: object=None
    planar_diagnostics: object=None

    def dissipation(self,velocity):
        result=super().dissipation(velocity)
        result['p9a_resistance_quadratic_W']=float(velocity@(self.planar_matrix@velocity))
        result['p9a_affine_power_W']=float(velocity@self.planar_rhs)
        result['full_resistance_quadratic_W']=float(velocity@(self.matrix@velocity))
        # Keep the old slip-dissipation diagnostic explicitly labelled, rather
        # than misrepresenting affine background power as passive dissipation.
        result['p65_slip_dissipation_W']=result.pop('total')
        return result


def augment_planar_system(system,shapes,mu,wall,gradient_at,gap_query):
    blocks=[];rhs=np.zeros_like(system.rhs);diagnostics=[]
    for k,i in enumerate(system.ids):
        shape=shapes[i]
        if type(shape) is not Sphere:raise TypeError('P9-A supports spherical MB only')
        if wall is None:
            blocks.append(np.zeros((6,6)));continue
        gap=gap_query(shape,wall)
        dr,db,diag=planar_wall_affine_block(shape.radius_m,mu,gap.gap_m,
            gap.normal_inward,gradient_at(shape.center_m),system.free[6*k:6*k+6])
        rim=rim_witness(shape,wall,gap)
        if rim is not None:
            dr[:]=0.;db[:]=0.
            diag.update(candidate_wall_weight=diag['wall_weight'],wall_weight=0.,
                open_rim_fallback=rim,applicability='OPEN_BOUNDARY_RIM_P65_FALLBACK')
        else:diag['applicability']='NEAREST_SINGLE_PLANAR_WALL'
        diag.update(particle_id=i,triangle_id=int(gap.wall_triangle_id),
            normal_role='AUTHORITATIVE_NEAREST_FINITE_WALL_FEATURE_NORMAL_LOCAL_PLANAR_APPROXIMATION')
        blocks.append(dr);rhs[6*k:6*k+6]=db;diagnostics.append(diag)
    delta=sparse.block_diag([sparse.csr_matrix(block) for block in blocks],format='csr')
    return PlanarResistanceSystem(system.ids,system.matrix+delta,system.rhs+rhs,
        system.self_diagonal,system.free,system.blocks,system.block_rows,delta,rhs,diagnostics)


class Particle9AStepper(SavedTrajectoryStepper):
    def __init__(self,*args,gradient_provider,**kwargs):
        super().__init__(*args,**kwargs)
        self.planar_statistics=dict(assembly_calls=0,wall_active_calls=0,
            maximum_abs_wall_normal_shear_force_N=0.,maximum_symmetry_error=0.)
        def transform(base,gap_query):
            def assembly(shapes,free,mu,wall=None,**options):
                system=augment_planar_system(base(shapes,free,mu,wall,**options),shapes,
                    mu,wall,gradient_provider,gap_query)
                stats=self.planar_statistics;stats['assembly_calls']+=1
                for d in system.planar_diagnostics:
                    stats['wall_active_calls']+=int(d['wall_weight']>0)
                    stats['maximum_abs_wall_normal_shear_force_N']=max(stats['maximum_abs_wall_normal_shear_force_N'],abs(d['wall_normal_shear_force_N']))
                    stats['maximum_symmetry_error']=max(stats['maximum_symmetry_error'],d['symmetry_error'])
                return system
            return assembly
        self.advance_cached,self.cache_info=cached_step_to(self,assembly_transform=transform)
        self.memoize=True
