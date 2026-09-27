"""Validation query candidates, separate from all physical eligibility decisions."""
from dataclasses import dataclass,asdict
import numpy as np
from .pair_broadphase import broadphase
from .pair_geometry import pair_gap
from .particle_shapes import Sphere
from .hydrodynamic_resistance import PhysicalNearField,evaluate_block


@dataclass(frozen=True)
class ValidationNeighborPolicy:
    center_cutoff_m:float
    skin_m:float
    reason:str
    role:str='VALIDATION_ONLY'

    def __post_init__(self):
        if not np.isfinite([self.center_cutoff_m,self.skin_m]).all() or self.center_cutoff_m<=0 or self.skin_m<0 or self.role!='VALIDATION_ONLY':
            raise ValueError('Only a finite positive validation neighbor policy is supported')

    @property
    def query_radius_m(self):return self.center_cutoff_m+self.skin_m

    def to_dict(self):return dict(asdict(self),cutoff_role='VALIDATION_NEIGHBOR_QUERY_ONLY',production_neighbor_cutoff_frozen=False,production_neighbor_skin_frozen=False)


def canonical_pairs(pairs):
    return sorted({tuple(sorted((int(i),int(j)))) for i,j in pairs if i!=j})


def standalone_candidates(shapes,policy):
    # Reuse the P4 AABB sweep; expand query bounds only, not physical shapes.
    raw=broadphase(shapes,surface_motion_bounds={i:policy.query_radius_m/2 for i in shapes})
    radius2=policy.query_radius_m**2
    return [(i,j) for i,j in raw if float(np.dot(shapes[i].center_m-shapes[j].center_m,shapes[i].center_m-shapes[j].center_m))<radius2]


def exact_physics_pairs(shapes,candidates,mu):
    records=[];selected=[]
    for i,j in canonical_pairs(candidates):
        g=pair_gap(shapes[i],shapes[j],i,j);active=g.state=='TOUCHING'
        reason='CONTACT_REGIME' if active else 'NONSPHERICAL_MODEL_NOT_FROZEN'
        if isinstance(shapes[i],Sphere) and isinstance(shapes[j],Sphere):
            spec=PhysicalNearField(i,j,g.gap_m,g.normal_j_to_i,g.roundoff_budget_m)
            _,_,d=evaluate_block(spec,shapes,mu);active=active or d['active'];reason=d['reason']
        if g.state=='PENETRATING':raise ValueError('Validation scene contains overlap')
        records.append(dict(pair=[i,j],gap_m=g.gap_m,shape_pair=[shapes[i].mode,shapes[j].mode],active=active,reason=reason))
        if active:selected.append((i,j))
    return selected,records
