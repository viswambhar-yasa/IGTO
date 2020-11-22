from geometry import Inputs
import numpy as np
import pytest


def test__controlpoints_coordinates_true():
    expected_output=np.array([[0., 0., 0., 1.], [2., 0., 0., 1.],[0., 2., 0., 1.],[2., 2., 0., 1.], [0., 0., 2., 1.], [2., 0., 2., 1.], [0., 2., 2., 1.], [2., 2., 2., 1.]])
    assert (np.array_equiv(Inputs(2,2,2,1,1,1).crtpts_coordinates(),expected_output)) is True


def test__knotvector_true():
    expected_output=np.array([0.,  0.,  0.,  0.,  0.5, 1.,  1.,  1.,  1. ])
    assert (np.array_equiv(Inputs().knot_vector(5,3),expected_output)) is True

