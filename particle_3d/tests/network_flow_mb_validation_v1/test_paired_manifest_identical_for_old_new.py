import csv,hashlib,json
import numpy as np
def test_paired_manifest_identical_for_old_new(report,events):
    p=report/'data/paired_cohort_manifest.csv';manifest=list(csv.DictReader(p.open()))
    assert len(manifest)==len(events)==30
    assert p.read_bytes()==(report/'server_bundle/data/paired_cohort_manifest.csv').read_bytes()
    for e in events:
        name=f"mb_{e['particle_id']:06d}"
        paths=[report/'outputs'/label/'trajectories'/name for label in ['OLD','NEW']]
        a,b=[np.load(p.with_suffix('.npz'))['samples'][0] for p in paths]
        assert a[:14].tobytes()==b[:14].tobytes()
        m,n=[json.loads(p.with_suffix('.json').read_text()) for p in paths]
        assert m['birth_metadata']==n['birth_metadata']==e
        assert m['integration_config']==n['integration_config']
        assert n['paired_identity']['cohort_manifest_sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
