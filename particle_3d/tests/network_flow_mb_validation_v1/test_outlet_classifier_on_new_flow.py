import json
def test_outlet_classifier_on_new_flow(report,newpaths,newenv):
    count=0
    for e,a in newpaths:
        m=json.loads((report/'outputs/NEW/trajectories'/f"mb_{e['particle_id']:06d}.json").read_text())
        if m['completed']:
            hit=newenv.classifier.first_event(a[-2,1:4],a[-1,1:4])
            assert hit is not None and hit.role==m['exit_outlet'];count+=1
    assert count>0
