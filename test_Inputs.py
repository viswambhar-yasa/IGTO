
from Preprocessing import Inputs
from numpy import array, array_equiv
import pytest

def test__controlpoints_coordinates_true():
    expected_output=array([[0., 0., 0., 1.], [2., 0., 0., 1.],[0., 2., 0., 1.],[2., 2., 0., 1.], [0., 0., 2., 1.], [2., 0., 2., 1.], [0., 2., 2., 1.], [2., 2., 2., 1.]])
    C=Inputs(2,2,2,2,2,2,1,1,1)
    actual_output=C.crtpts_coordinates()
    assert (array_equiv(actual_output[:,:-1],expected_output)) is True


def test__knotvector_true():
    '''
    Values obtained from NURBS book
    '''
    expected_output=array([0.,  0.,  0.,  0.,  0.5, 1.,  1.,  1.,  1. ])
    assert (array_equiv(Inputs().knot_vector(5,3),expected_output)) is True

def test__knotconnect_span_true():
    '''
    Values obtained from NURBS book
    '''
    C=Inputs(1,1,1,5,5,5,1,1,1)
    XI_SPAN,XI_KNOTCONNECTIVITY,XI_UNIKNOTS,nU=C.xi_knotspan()
    expected_output=array([[0,0.25],[0.25, 0.5 ],[0.5,0.75],[0.75, 1]])
    assert array_equiv(XI_SPAN,expected_output) is True