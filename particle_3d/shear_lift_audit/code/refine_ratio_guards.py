"""Apply operator-specific roundoff budgets to already audited states; no FEM resampling.

The lubrication operator can amplify normal-velocity roundoff. Its denominator
budget must therefore use its own coefficient, not only the self-resistance.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
from audit_math import candidate,force_ratio
from run_audit import aggregate,dump

def main(config):
    cfg=json.loads(Path(config).read_text());out=Path(cfg['output']);receipts=json.loads((out/'trajectory_audit_receipts.json').read_text())
    changes=0;rows=0;physical_changed=False
    for r in receipts:
        if not r['count']:continue
        p=Path(r['file'])
        with np.load(p) as z:d={k:z[k] for k in z.files}
        physical_before=hashlib.sha256(d['CANDIDATE_SAFFMAN_MAGNITUDE'].tobytes()+d['drag_N'].tobytes()+d['lubrication_N'].tobytes()).hexdigest()
        diag=6*np.pi*.00345312*d['radius_m'];condition=1+d['lubrication_coefficient_kg_s']/diag
        velocity_zero=512*6*np.finfo(float).eps*condition*np.maximum(np.linalg.norm(d['local_u_m_s'],axis=1),np.linalg.norm(d['particle_v_m_s'],axis=1))
        d['velocity_numerical_zero_m_s']=velocity_zero;d['force_numerical_zero_N']=diag*velocity_zero
        d['self_scaled_resistance_condition']=condition
        d['candidate_direction_only'][d['slip_speed_m_s']<=velocity_zero]=np.nan
        lift_zero=candidate(d['radius_m'],velocity_zero,d['shear_s_inv'],.00345312,1056.)
        lub_zero=np.maximum(d['force_numerical_zero_N'],d['lubrication_coefficient_kg_s']*velocity_zero)
        ld,ls=force_ratio(d['CANDIDATE_SAFFMAN_MAGNITUDE'],d['drag_N'],d['force_numerical_zero_N'],lift_zero)
        ll,ss=force_ratio(d['CANDIDATE_SAFFMAN_MAGNITUDE'],d['lubrication_N'],lub_zero,lift_zero)
        ll[d['lubrication_coefficient_kg_s']==0]=np.nan;ss[d['lubrication_coefficient_kg_s']==0]=3
        changes+=int(np.sum(ss!=d['lift_lubrication_status']));rows+=r['count']
        d.update(lift_drag_ratio=ld,lift_drag_status=ls,lift_lubrication_ratio=ll,lift_lubrication_status=ss,
                 lift_numerical_zero_N=lift_zero,lubrication_numerical_zero_N=lub_zero)
        physical_after=hashlib.sha256(d['CANDIDATE_SAFFMAN_MAGNITUDE'].tobytes()+d['drag_N'].tobytes()+d['lubrication_N'].tobytes()).hexdigest()
        assert physical_before==physical_after
        np.savez_compressed(p,**d)
    aggregate(receipts,cfg)
    dump(out/'ratio_guard_refinement.json',dict(status='PASS',state_count=rows,changed_lubrication_ratio_statuses=changes,
         reason='OPERATOR_SPECIFIC_PROPAGATION_OF_VELOCITY_ROUNDOFF',force_arrays_changed=False,FEM_resampled=False,trajectory_recomputed=False,
         denominator_budget='max(Stokes_drag_budget, lubrication_coefficient * velocity_roundoff_budget)',
         numerator_budget='Saffman_prefactor * velocity_roundoff_budget'))
    print('RATIO_GUARD_PASS',rows,changes)

if __name__=='__main__':main(sys.argv[1])
