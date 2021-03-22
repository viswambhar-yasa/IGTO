#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#------------------------------------#
#A python test file to check input.py 
# command to run all test cases
# pytest test_Inputs.py
# -----------------------------------# 

from Preprocessing import Inputs
from numpy import array, array_equiv
import pytest

def test__controlpoints_coordinates_true():
    '''
    UNIT TESTING
    Aim: The co-ordinates of the control points are generated for crtpts_coordinates method in Input class

    Expected result : Values of a cube with equal length,width and height(a total of 8 co-ordinate points have to be generated)

    Test command : pytest test_Inputs.py::test__controlpoints_coordinates_true

    Remarks : test case passed successfully
    '''
    expected_output=array([[0., 0., 0., 1.], [2., 0., 0., 1.],[0., 2., 0., 1.],[2., 2., 0., 1.], [0., 0., 2., 1.], [2., 0., 2., 1.], [0., 2., 2., 1.], [2., 2., 2., 1.]])
    C=Inputs(2,2,2,2,2,2,1,1,1)
    actual_output=C.crtpts_coordinates()
    assert (array_equiv(actual_output[:,:-1],expected_output)) is True


def test__knotvector_true():
    '''
    UNIT TESTING
    Aim: The knot vector is generated for number of control points and degree of the basis from knot_vector method from input class

    Expected result : For 5 control points and a cubic basis. The knot vector is given as [0.,  0.,  0.,  0.,  0.5, 1.,  1.,  1.,  1. ]
                        The results are obtained from NURBS book
    
    Test command : pytest test_Inputs.py::test__knotvector_true

    Remarks : test case passed successfully
    '''
    expected_output=array([0.,  0.,  0.,  0.,  0.5, 1.,  1.,  1.,  1. ])
    assert (array_equiv(Inputs().knot_vector(5,3),expected_output)) is True

def test__knotconnect_span_true():
    '''
    UNIT TESTING
    Aim: The length of each element is represented through knot span. Knot span is defined as a method from input class

    Expected result : For given inputs values and 5 number of element along x,y and z direction. The knot span along x direction is given as [[0,0.25],[0.25, 0.5 ],[0.5,0.75],[0.75, 1]]
                        The results are obtained from NURBS book
    
    Test command : pytest test_Inputs.py::test__knotconnect_span_true
    
    Remarks : test case passed successfully
    '''
    C=Inputs(1,1,1,5,5,5,1,1,1)
    XI_SPAN,XI_KNOTCONNECTIVITY,XI_UNIKNOTS,nU=C.xi_knotspan()
    expected_output=array([[0,0.25],[0.25, 0.5 ],[0.5,0.75],[0.75, 1]])
    assert array_equiv(XI_SPAN,expected_output) is True