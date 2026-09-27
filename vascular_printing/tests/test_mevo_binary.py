"""Retained binary membership and UNKNOWN accounting regressions."""
import numpy as np
import pytest
from vascular_processing.label_transfer_qc import collapse_labels_to_mevo,evaluate_mevo_binary,evaluate_mevo_boundary

def test_native_membership_merge_preserves_unknown_and_internal_errors():
    truth=np.array([1,2,3,2,3]);prediction=np.array([1,3,2,0,0])
    np.testing.assert_array_equal(collapse_labels_to_mevo(truth,truth=True),[1,2,2,2,2])
    np.testing.assert_array_equal(collapse_labels_to_mevo(prediction),[1,2,2,0,0])
    result=evaluate_mevo_binary(truth,prediction)
    assert result['within_mevo']['now_correct_binary_count']==2
    assert result['MeVO_precision']==1
    assert result['MeVO_recall']==.5
    assert result['MeVO_f1']==pytest.approx(2/3)
    np.testing.assert_array_equal(prediction,[1,3,2,0,0])
    with pytest.raises(ValueError):collapse_labels_to_mevo([0,1],truth=True)


def test_false_inclusion_and_both_exclusions_use_full_truth_denominators():
    result=evaluate_mevo_binary([1,1,1,2,2,3,3],[1,2,0,3,1,2,0])
    assert result['confusion_matrix']==[[1,1,1],[1,2,1]]
    assert result['M1_precision']==.5
    assert result['M1_recall']==pytest.approx(1/3)
    assert result['M1_f1']==pytest.approx(.4)
    assert result['MeVO_f1']==pytest.approx(4/7)
    assert result['false_inclusion']==pytest.approx(1/3)
    assert result['false_exclusion_to_M1']==.25
    assert result['false_exclusion_to_UNKNOWN']==.25
    assert result['false_exclusion_total']==.5
    assert result['UNKNOWN_fraction']==pytest.approx(2/7)


def test_internal_m2_m3_swap_does_not_move_proximal_boundary():
    points=np.c_[np.arange(12),np.zeros((12,2))]
    truth=np.repeat([1,2,3],4);prediction=np.repeat([1,3,2],4)
    assert evaluate_mevo_boundary(points,truth,prediction)==0
    assert evaluate_mevo_boundary(points,truth,np.zeros(12,dtype=int)) is None

