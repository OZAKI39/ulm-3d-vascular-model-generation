"""Stable NN sampling and source protection checks."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from vascular_processing.topbrain_semantic_points import Geometry,SemanticPoints,downsample,trunk_direction
from vascular_processing.similarity_registration import transform_points
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

