"""Independent objectives and explicit constraints on actual material components."""
from common import *

BINS=[(1,2),(3,13),(14,33),(34,113),(114,200),(201,None)]
SPEC=dict(
 name='EXPLORATORY_MULTIFRAGMENT_TARGET',target_nodes=[20,200],min_target_count=5,
 min_total_detached_fraction_initial=0.,max_largest_fraction_detached=.50,
 max_fine_debris_fraction_detached=.20,post_unload_connectivity_required=True,
 common_cohort_time_s=.0034,observation_end_s=END,
 size_bins=BINS,
 meaning='Exploratory engineering constraints at this discrete resolution, not clinical size classes or validated biological criteria. V2 removes the old extra 10% mass floor; its specified count, dominance and fines constraints are used.',
 optimization='Return separate objectives, constraints and per-fragment evidence; no weighted sum hides a failed constraint.',
)

def shape(x,dx=DX):
    count=len(x);center=x.mean(0);q=x-center
    covariance=q.T@q/count;values,axes=np.linalg.eigh(covariance);order=np.argsort(values)[::-1];values=np.maximum(values[order],0);axes=axes[:,order];rms=np.sqrt(values)
    tol=1e-8*dx;rank=int(np.sum(rms>tol))
    ratios=[float(rms[0]/v) if v>tol else None for v in rms]
    label='UNRESOLVED_FINE_DEBRIS' if count<=2 else 'CHAIN_LIKE' if rank<2 else 'SHEET_LIKE' if rank<3 or rms[2]<.15*rms[0] else 'THREE_DIMENSIONAL_POINT_CLOUD'
    return dict(actual_span_m=np.ptp(x,axis=0),radius_gyration_m=float(np.sqrt(values.sum())),principal_axes_columns=axes,principal_rms_lengths_m=rms,length_to_principal_ratios=ratios,point_cloud_rank=rank,shape_class=label,definition='Measured point coordinates, no smoothed shell, no invented thickness; degenerate ratios are undefined')

def metrics(sizes,n):
    sizes=np.asarray(sizes,int);detached=int(sizes.sum());largest=int(sizes.max(initial=0));target=(sizes>=20)&(sizes<=200);fine=int(sizes[sizes<=2].sum());fines=fine/detached if detached else None;dominance=largest/detached if detached else None
    criteria=dict(target_count=int(target.sum())>=SPEC['min_target_count'],detached_mass=detached/n>=SPEC['min_total_detached_fraction_initial'],no_dominant_fragment=dominance is not None and dominance<=SPEC['max_largest_fraction_detached'],fine_debris=fines is not None and fines<=SPEC['max_fine_debris_fraction_detached'])
    out=dict(fragment_count_target_band=int(target.sum()),mass_fraction_in_target_band=float(sizes[target].sum()/n),largest_free_fragment_fraction_initial=largest/n,largest_free_fragment_fraction_detached=dominance,fine_debris_mass_fraction_detached=fines,fine_debris_mass_fraction_initial=fine/n,total_detached_mass_fraction=detached/n,attached_body_fraction=1-detached/n,N_ge2=int(np.sum(sizes>=2)),N_ge14=int(np.sum(sizes>=14)),N_ge20=int(np.sum(sizes>=20)),N_ge34=int(np.sum(sizes>=34)),N_ge114=int(np.sum(sizes>=114)),N_ge200um=int(np.sum((6*sizes*V/np.pi)**(1/3)>=200e-6)),N_ge300um=int(np.sum((6*sizes*V/np.pi)**(1/3)>=300e-6)),instant_target=all(criteria.values()),independent_criteria=criteria)
    out.update(n_fragments_20_200=int(target.sum()),n_fragments_34_200=int(np.sum((sizes>=34)&(sizes<=200))),n_fragments_ge_114=int(np.sum(sizes>=114)),mass_fraction_target_initial=float(sizes[target].sum()/n),mass_fraction_target_detached=float(sizes[target].sum()/detached) if detached else None)
    return out

def scalar_metrics(sizes,n):
    return {key:value for key,value in metrics(sizes,n).items() if key!='independent_criteria'}

def selection_key(sizes,n):
    q=metrics(sizes,n)
    return (q['fragment_count_target_band'],q['mass_fraction_in_target_band'],q['N_ge14'],q['total_detached_mass_fraction'])

def freeze():
    path=O/'MULTIFRAGMENT_OBJECTIVE.json'
    assert not path.exists(),'Do not change a frozen objective based on outcomes'
    dump(path,SPEC)

if __name__=='__main__':freeze()
