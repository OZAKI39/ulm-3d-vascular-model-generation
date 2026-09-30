import json
from pathlib import Path
import numpy as np
from pd_clot.geometry import make_cloud
from pd_clot.fragment_topology import FragmentTracker
from pd_clot.fragment_mechanics import prepare_supported_shape
from pd_clot.regularization.diagnostics import classify, ClassifiedClearance, localization


def test_components_are_classified_without_merging_or_deletion():
    c=json.loads((Path(__file__).parents[1]/'configs/streaming_regularized_demo_coarse.json').read_text())
    g=make_cloud(c['clot']);b=np.ones(len(g.pairs));tracker=FragmentTracker(g,.001)
    tracker.update(b,g.X,np.zeros_like(g.X),0)
    # Graph diagnostic fixture, not a physical run or a fracture driver.
    victim=int(np.flatnonzero(~g.fixed)[0]);b[np.any(g.pairs==victim,axis=1)]=0
    attached,comps,_=tracker.update(b,g.X,np.zeros_like(g.X),1)
    ranks=prepare_supported_shape(g,b,c['safety'])[3]
    category,rows,s=classify(g,b,ranks,attached,tracker.labels,comps,c['fragment_resolution'])
    assert category[victim]==3 and s['P_singleton']==1 and s['P_lowrank']==1
    assert s['fragment_resolution_status']=='FRAGMENT_RESOLUTION_INADEQUATE'
    assert sum(r['particle_count'] for r in rows)==len(g.X)
    clear=ClassifiedClearance(len(g.X));mask=np.zeros(len(g.X),bool);mask[victim]=True
    first=clear.update(mask,category,g.volume);second=clear.update(mask,category,g.volume)
    assert first==second and first['resolved_fragment_clearance_fraction']==0
    assert np.isclose(first['under_resolved_debris_clearance_fraction'],1/len(g.X))
    loc=localization(g,np.zeros(len(g.X)),b,attached,np.array(c['streaming']['bubble_center_m']),.000525)
    assert loc['fraction_of_damage_inside_influence_region'] is None
