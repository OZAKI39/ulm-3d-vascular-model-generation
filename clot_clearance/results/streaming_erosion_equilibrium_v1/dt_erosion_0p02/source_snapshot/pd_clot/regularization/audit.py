"""Read-only scientific identity and state audits for a completed candidate."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def audit_run(path):
    p=Path(path);z=np.load(p/'states.npz');h=json.loads((p/'history.json').read_text());c=json.loads((p/'CONFIG.json').read_text())
    summary=json.loads((p/'SUMMARY.json').read_text());pairs=z['pairs'];volume=z['volume'];fixed=z['fixed'];N=z['cycles']
    checks={'finite_saved_arrays':all(np.isfinite(z[k]).all() for k in z.files),
            'full_exposure':N[-1]==c['damage']['DeltaN']*c['simulation']['number_of_macro_steps'],
            'damage_monotonic':bool(np.all(np.diff(z['bond_D'],axis=0)>=0)),
            'D_break_one':c['damage']['D_break']==1 and summary['D_break']==1,
            'integrity_matches_damage':np.array_equal(z['integrity'],1-z['bond_D']),
            'active_bonds_match_complete_damage':np.array_equal(z['bond_active'],z['bond_D']<1),
            'all_particles_retained':z['x'].shape[1]==len(z['X']),
            'fixed_positions_exact':bool(np.all(z['x'][:,fixed]==z['X'][fixed])),
            'positive_jacobian':summary['minimum_J_or_intrinsic_stretch']>c['safety']['minimum_J'],
            'source_is_upstream':c['streaming']['bubble_center_m'][0]<c['clot']['origin_m'][0]}
    graph=True;energy=True;classification=True
    for k,row in enumerate(h):
        q=pairs[z['bond_active'][k]];i=q.ravel();j=q[:,::-1].ravel()
        _,labels=connected_components(coo_matrix((np.ones(len(i)),(i,j)),shape=(len(volume),len(volume))),directed=False)
        attached=np.isin(labels,np.unique(labels[fixed]));graph &= np.array_equal(attached,z['attached'][k])
        energy &= bool(np.isclose(z['bond_dissipated_energy_J'][k].sum(),row['damage_dissipation_estimate_J'],rtol=1e-9,atol=1e-20))
        energy &= abs(row['relative_numerical_energy_residual'])<=c['energy_audit']['maximum_relative_residual'] or abs(row['numerical_energy_residual_J'])<=c['energy_audit']['minimum_absolute_tolerance_J']
        freevol=volume[~attached].sum();classes=z['resolution_class'][k]
        classification &= bool(np.isclose(sum(volume[classes==v].sum() for v in [1,2,3]),freevol,rtol=1e-12,atol=1e-30))
    checks.update(graph_reconstructed=graph,energy_ledger_closes=energy,detached_categories_conserve_volume=classification)
    before=json.loads((p/'IDENTITY.json').read_text())['source_sha256']
    checks['run_source_snapshot_identity']=all(hashlib.sha256((p/'source_snapshot'/name).read_bytes()).hexdigest()==sha for name,sha in before.items())
    checks={k:bool(v) for k,v in checks.items()}
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
            'first_failure':summary['events']['first_bond_failure'],'first_detachment':summary['events']['first_detachment'],
            'maximum_energy_residual':max(abs(r['relative_numerical_energy_residual']) for r in h),
            'note':'Software/ledger verification only; not physical validation or mesh convergence.'}
    (p/'AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    if result['status']!='PASS':raise AssertionError(result)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args();print(json.dumps(audit_run(a.run),indent=2))


if __name__=='__main__':main()
