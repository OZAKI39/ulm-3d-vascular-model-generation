"""Re-evaluate saved V2 decomposition arrays without solving or changing gates."""
import numpy as np
from ..paths import ROOT, read_json, write_json, fresh_directory, sha256
from .field_error_decomposition import decompose, classification

def review(output, root=ROOT):
    output=fresh_directory(output,root)
    path=root/'diagnostics/reference_pipe/gmsh_reference/error_diagnosis/DECOMPOSITION_SAMPLES.npz'
    with np.load(path) as z:
        groups=[]
        for group in np.unique(z['groups']):
            m=z['groups']==group
            d=decompose(z['u_CFD'][m],z['u_P1'][m],z['u_exact'][m])
            groups.append(dict(group=int(group),classification=classification(d),velocity=d))
    original=read_json(root/'diagnostics/reference_pipe/V2_FIELD_ERROR_DECOMPOSITION.json')
    assert len(groups)==len(original['groups'])
    for actual,expected in zip(groups,original['groups']):
        for k in ['total_relative_L2','space_relative_L2','discrete_relative_L2']:
            np.testing.assert_allclose(actual['velocity'][k],expected['velocity'][k],rtol=1e-12,atol=1e-14)
        assert actual['classification']==expected['classification']
    result=dict(status='PASS_REPRODUCED_SAVED_DIAGNOSTIC',scientific_state='MIXED_FIELD_ERROR',
                input_sha256=sha256(path),groups=groups,new_CFD_runs=0,new_particle_runs=0)
    write_json(output/'REVIEW.json',result)
    return result
