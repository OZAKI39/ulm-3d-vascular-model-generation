"""Behavioral checks for proper registration, soft transport and data protection."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from vascular_processing.topbrain_semantic_points import Geometry,SemanticPoints,downsample,trunk_direction
from vascular_processing.similarity_registration import register,transform_points,proper_similarity
from vascular_processing.partial_fgw_transfer import transfer,transport_probabilities
from vascular_processing.semantic_ensemble import ensemble
from vascular_processing.label_transfer_qc import scores

ROOT=Path(__file__).resolve().parents[1]


def fixture():
    t=np.linspace(0,12,24)
    pieces=[np.c_[t,.1*np.sin(t),.03*t]];classes=[np.ones(len(t),int)]
    for sign in [-1,1]:
        u=np.linspace(.3,1,24)
        start=pieces[0][-1]
        m2=start+u[:,None]*np.array([13,sign*9,sign*3])
        pieces.append(m2);classes.append(np.full(len(u),2))
        for daughter in [-1,1]:
            v=np.linspace(.05,1,24)
            m3=m2[-1]+v[:,None]*np.array([8+3*daughter,sign*(6+2*daughter),daughter*9+sign*4])
            m3[:,2]+=1.5*np.sin(v*4)
            pieces.append(m3);classes.append(np.full(len(v),3))
    points=np.vstack(pieces);labels=np.concatenate(classes)
    g=Geometry(points,np.ones(len(points))/len(points),points[0],trunk_direction(points,points[0]),'donor','L')
    donor=SemanticPoints(g,labels,points,labels,np.arange(len(points)),np.arange(len(points)),{})
    matrix=np.eye(4);matrix[:3,:3]=1.4*Rotation.from_rotvec([.2,.3,.7]).as_matrix();matrix[:3,3]=[12,-5,7]
    target_points=transform_points(points,matrix)
    target=Geometry(target_points,g.mass.copy(),target_points[0],matrix[:3,:3]@g.direction/1.4,'target','L')
    return donor,target,matrix


def test_similarity_transform_recovery_and_synthetic_semantics():
    donor,target,matrix=fixture();registration=register(donor.geometry,target)
    assert registration['status']=='VALID'
    # ICP is coarse alignment of unpaired samples, not exact paired Procrustes.
    error=transform_points(donor.geometry.points,registration['transform'])-target.points
    size=np.percentile(np.linalg.norm(target.points-target.root,axis=1),90)
    assert np.sqrt(np.mean(np.sum(error**2,axis=1)))/size < .02
    scale=proper_similarity(registration['transform'])
    assert scale==pytest.approx(1.4,rel=.03)
    rotation_error=Rotation.from_matrix(registration['transform'][:3,:3]/scale@(matrix[:3,:3]/1.4).T).magnitude()
    assert np.rad2deg(rotation_error)<2
    result=transfer(donor,target,registration)
    assert scores(donor.labels,result['predicted_label'])['macro_f1']>.80
    assert result['transported_mass']==pytest.approx(.8)


def test_reflection_and_side_mismatch_rejected():
    matrix=np.diag([-1.,1,1,1])
    with pytest.raises(ValueError,match='Reflection'):proper_similarity(matrix)
    with pytest.raises(ValueError,match='Nonfinite'):proper_similarity(np.full((4,4),np.nan))
    with pytest.raises(ValueError,match='collapsed'):proper_similarity(np.diag([0.,0,0,1]))
    donor,target,_=fixture();target.side='R'
    with pytest.raises(ValueError,match='SIDE_MISMATCH'):register(donor.geometry,target)


def test_transport_mass_soft_semantics_and_unknown():
    t=np.array([[.1,.1,0,.00001],[0,.3,0,0],[0,0,0,0]])
    r=transport_probabilities(t,np.array([1,2,3]))
    np.testing.assert_allclose(r['probabilities'][1],[.25,.75,0])
    assert r['predicted_label'].tolist()==[1,2,0,0]


def test_equal_ensemble_counts_supported_donors():
    a=transport_probabilities(np.array([[.2,0],[0,.3],[0,0]]),[1,2,3])
    b=transport_probabilities(np.array([[0,0],[.2,0],[0,.3]]),[1,2,3])
    c=ensemble([a,b])
    np.testing.assert_allclose(c['probabilities'],[[.5,.5,0],[0,.5,.5]])
    assert c['valid_donor_count'].tolist()==[2,2]
    empty=transport_probabilities(np.zeros((3,2)),[1,2,3])
    assert ensemble([a,empty,empty])['predicted_label'].tolist()==[0,0]


def test_extra_distal_geometry_can_remain_unknown():
    donor,target,matrix=fixture()
    extra=target.points[-1]+np.c_[np.linspace(30,50,32),np.linspace(20,40,32),np.linspace(-20,-40,32)]
    target.points=np.vstack([target.points,extra]);target.mass=np.ones(len(target.points))/len(target.points)
    registration=dict(status='VALID',transform=matrix)
    result=transfer(donor,target,registration)
    assert np.mean(result['predicted_label'][-32:]==0)>=.5
    assert scores(donor.labels,result['predicted_label'][:len(donor.labels)])['macro_f1']>.75


def test_downsampling_preserves_actual_samples_labels_and_mass():
    donor,_,_=fixture();before=donor.native_points.copy()
    ids,assignment,mass=downsample(donor.native_points,35)
    assert len(ids)<=35 and len(assignment)==len(before)
    assert mass.sum()==pytest.approx(1.)
    assert np.all(assignment[ids]==np.arange(len(ids)))
    np.testing.assert_array_equal(donor.native_points,before)
    assert set(donor.native_labels[ids])=={1,2,3}


def test_source_hashes_and_s1_2_protected():
    baseline=json.loads((ROOT/'outputs/topbrain_brava_transfer/protection_before.json').read_text())
    for name,expected in baseline.items():
        if name.endswith('.nii.gz'):continue  # Full NIfTI scan is independently recorded by the final audit.
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    assert baseline['s1-2_swc_roi_generate_human.py']=='3361cea44a97e8bfa907ff1504b0f977e329d35eb59cb9f57cc1b3a2a55cd729'


def test_absolute_f1_pass_does_not_override_baseline_stop_rule():
    from vascular_processing.label_transfer_qc import pilot_gate
    assert not pilot_gate(dict(macro_f1=.77,M2_f1=.84,M3_f1=.70),dict(macro_f1=.775))
    assert pilot_gate(dict(macro_f1=.78,M2_f1=.84,M3_f1=.70),dict(macro_f1=.775))
