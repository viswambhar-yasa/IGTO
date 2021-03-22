#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#---------------------------------------#
#A python test file to check geomeetry.py 
# command to run all test cases
# pytest test_geometry.py
# --------------------------------------# 
from Preprocessing import Inputs
from geometry import knot_index,bspline_basis,derbspline_basis,trilinear_der,controlpointassembly
from numpy import array,array_equiv,sum,equal,round,ones,float,zeros,all
import pytest

def test__knot_index_true():
    '''
    UNIT TESTING
    Aim: The knot index has to be calculated for a given knot vector, the position or index of the value nearest ot given knot should be returned. 

    Expected result : For given knotvector=[0,0,0,0.14285714,0.28571429,0.42857143,0.57142857,0.71428571,0.85714286,1,1,1]
                            degree=2
                            U=0.3
                        The results are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__knot_index_true
    
    Remarks : test case passed successfully
    '''
    knotvector=[0,0,0,0.14285714,0.28571429,0.42857143,0.57142857,0.71428571,0.85714286,1,1,1]
    degree=2
    U=0.3
    assert (knot_index(degree,U,knotvector)==4) is True

def test__Bspline_basis_equal_true():
    '''
    UNIT TESTING
    Aim: The Bspline basis are calculated for given knot vector and degree and knot value

    Expected result : For given knotvector=[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
                            degree=2
                            U=5/2
                            BASIS : [0.125,0.75,0.125]
                        The results are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__Bspline_basis_equal_true
    
    Remarks : test case passed successfully
    '''
    knotvector=[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
    degree=2
    U=5/2
    knotindex=knot_index(degree,U,knotvector)
    output=bspline_basis(knotindex,degree,U,knotvector)
    assert (array_equiv(output,[0.125,0.75,0.125])) is True

def test__Bspline_basis_sum_equal_to_one_true():
    '''
    UNIT TESTING
    Aim: The Bspline basis property needs to be satisfied which is the sum of all basis function is unity.

    Expected result : For given  knotvector=[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
                                    degree=2
                                    U=5/3
                        The knotvector,degree and knot value are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__Bspline_basis_sum_equal_to_one_true

    Remarks : test case passed successfully
    '''
    knotvector=[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
    degree=2
    U=5/3
    knotindex=knot_index(degree,U,knotvector)
    output=sum(bspline_basis(knotindex,degree,U,knotvector))
    assert (int(output)==1) is True

def test__Bspline_basis_all_values_positive_true():
    '''
    UNIT TESTING
    Aim: The Bspline basis property needs to be satisfied which is the sum of all basis function is unity.

    Expected result : For given  knotvector=[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
                                    degree=2
                                    U=5/3
                                    Sum of all basis should be 1
                        The knotvector,degree and knot value are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__Bspline_basis_all_values_positive_true

    Remarks : test case passed successfully
    '''
    knotvector=[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
    degree=2
    U=5/3
    knotindex=knot_index(degree,U,knotvector)
    output=all(bspline_basis(knotindex,degree,U,knotvector)>=0)
    o=0
    if output:
        o=1
    assert (o==1) is True


def test__derbspline_basis_equal_true():
    '''
    UNIT TESTING
    Aim: The derivatives Bspline basis are obtain for the given knotvecor,degree and knot value.

    Expected result : For given   knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
                                    degree=3
                                    U=0.6
                                Derivative of basis : [-0.72, -2.24,  2.48,  0.48]
                        The knotvector,degree and knot value are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__derbspline_basis_equal_true

    Remarks : test case passed successfully
    '''
    knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
    degree=3
    U=0.6
    knotindex=knot_index(degree,U,knotvector)
    output=round(derbspline_basis(knotindex,degree,U,knotvector),2 )
    expected_output=array([-0.72, -2.24,  2.48,  0.48])
    assert (array_equiv(output,expected_output)) is True

def test__derbspline_basis_sum_equal_zero_true():
    '''
    UNIT TESTING
    Aim: The Derivatives of Bspline basis property needs to be satisfied which is the sum of all basis function is Zero.

    Expected result : For given   knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
                                    degree=3
                                    U=0.6
                                    Sum of all derivatives should be 0
                        The knotvector,degree and knot value are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__derbspline_basis_sum_equal_zero_true

    Remarks : test case passed successfully
    '''
    knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
    degree=3
    U=0.6
    knotindex=knot_index(degree,U,knotvector)
    output=round(sum(derbspline_basis(knotindex,degree,U,knotvector)),8)
    assert (float(output)==0.0) is True



def test__trilinear_der_Basis_sum_equal_to_one_true():
    '''
    UNIT TESTING
    Aim: Trivariant function are calculated for given knot vectors and it should satisfy the property of bspline basis i.e sum of basis is unity

    Expected result : For given   knotvector=[0,0,0.5,1,1]
                                    degree=1
                                    U=0.6,0.5,0.3
                                    Sum of basis  should be 1
                        The knotvectors,degree and knot values are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__trilinear_der_Basis_sum_equal_to_one_true

    Remarks : test case passed successfully
    '''
    WEIGHTS=ones(8)
    Xi_degree=Eta_degree=Neta_degree=1
    Xi_knotvector=Eta_knotvector=Eta_knotvector=[0,0,0.5,1,1]
    DR_DX,DR_DY,DR_DZ,R =trilinear_der(0.6,0.5,0.3,WEIGHTS,Xi_degree,Xi_knotvector,Eta_degree,Eta_knotvector,Neta_degree,Eta_knotvector)
    output=sum(R)
    assert (float(output)==1.0) is True


def test__trilinear_der_Basis_less_than_zero_false():
    '''
    UNIT TESTING
    Aim: Trivariant function are calculated for given knot vectors and it should satisfy the property of bspline basis i.e basis should not be negative

    Expected result : For given   knotvector=[0,0,0,0,0,0.16666667,0.33333333,0.5,0.6666666,0.83333333,1,1,1,1,1]
                                    degree=4
                                    U=0.6,0.5,0.3
                                    basis should not be negative
                        The knotvectors,degree and knot values are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__trilinear_der_Basis_less_than_zero_false

    Remarks : test case passed successfully
    '''
    
    Xi_degree=Eta_degree=Neta_degree=4
    WEIGHTS=ones((Xi_degree+1)**3)
    Xi_knotvector=Eta_knotvector=Eta_knotvector=[0,0,0,0,0,0.16666667,0.33333333,0.5,0.6666666,0.83333333,1,1,1,1,1]
    DR_DX,DR_DY,DR_DZ,R =trilinear_der(0.6,0.5,0.3,WEIGHTS,Xi_degree,Xi_knotvector,Eta_degree,Eta_knotvector,Neta_degree,Eta_knotvector)
    output=all(R<=0)
    o=0
    if output==True:    
        o=1
    assert (o==1) is False

def test__trilinear_der_XI_sum_equal_to_zero_true():
    '''
    UNIT TESTING
    Aim: Derivative of Trivariant function are calculated for given knot vectors and it should satisfy the property of derivative of bspline basis i.e sum of derivative basis should not be zero

    Expected result : For given   knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
                                    degree=2
                                    U=1/4,2/3,3/5
                                    sum of derivative of trivariant function w.r.t to xi should be zero
                        The knotvectors,degree and knot values are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__trilinear_der_XI_sum_equal_to_zero_true

    Remarks : test case passed successfully
    '''
    Xi_degree=Eta_degree=Neta_degree=2
    WEIGHTS=ones((Xi_degree+1)**3)
    Xi_knotvector=Eta_knotvector=Eta_knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
    DR_DXI,DR_DETA,DR_DNETA,R =trilinear_der(1/4,2/3,3/5,WEIGHTS,Xi_degree,Xi_knotvector,Eta_degree,Eta_knotvector,Neta_degree,Eta_knotvector)
    output=float(round(sum(DR_DXI),10))
    assert (output==0.0) is True


def test__trilinear_der_ETA_sum_equal_to_zero_true():
    '''
    UNIT TESTING
    Aim: Derivative of Trivariant function are calculated for given knot vectors and it should satisfy the property of derivative of bspline basis i.e sum of derivative basis should not be zero

    Expected result : For given   knotvector=[0,0,0,0.16666667,0.33333333, 0.5,0.66666667, 0.83333333,1,1,1]
                                    degree=2
                                    U=0.5,0.8,0.2
                                    sum of derivative of trivariant function w.r.t to eta should be zero
                        The knotvectors,degree and knot values are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__trilinear_der_ETA_sum_equal_to_zero_true

    Remarks : test case passed successfully
    '''
    Xi_degree=Eta_degree=Neta_degree=2
    WEIGHTS=ones((Xi_degree+1)**3)
    Xi_knotvector=Eta_knotvector=Eta_knotvector=[0,0,0,0.16666667,0.33333333, 0.5,0.66666667, 0.83333333,1,1,1]
    DR_DXI,DR_DETA,DR_DNETA,R =trilinear_der(0.5,0.8,0.2,WEIGHTS,Xi_degree,Xi_knotvector,Eta_degree,Eta_knotvector,Neta_degree,Eta_knotvector)
    output=float(round(sum(DR_DETA),10))
    assert(output==0) is True

def test__trilinear_der_NETA_sum_equal_to_zero_true():
    '''
    UNIT TESTING
    Aim: Derivative of Trivariant function are calculated for given knot vectors and it should satisfy the property of derivative of bspline basis i.e sum of derivative basis should not be zero

    Expected result : For given   knotvector=[0,0,0,0,0.33333333, 0.66666667,1,1,1,1]
                                    degree=3
                                    U=0.5,0.8,0.2
                                    sum of derivative of trivariant function w.r.t to neta should be zero
                        The knotvectors,degree and knot values are obtained from NURBS book
    
    Test command : pytest test_geometry.py::test__trilinear_der_NETA_sum_equal_to_zero_true

    Remarks : test case passed successfully
    '''
    Xi_degree=Eta_degree=Neta_degree=3
    WEIGHTS=ones((Xi_degree+1)**3)
    Xi_knotvector=Eta_knotvector=Eta_knotvector=[0,0,0,0,0.33333333, 0.66666667,1,1,1,1]
    DR_DXI,DR_DETA,DR_DNETA,R =trilinear_der(0.5,0.8,0.2,WEIGHTS,Xi_degree,Xi_knotvector,Eta_degree,Eta_knotvector,Neta_degree,Eta_knotvector)
    output=float(round(sum(DR_DNETA),10))
    assert(output==0) is True


def test__single_element_Assembly():
    '''
    UNIT TESTING
    Aim: Control point assembly is tested for a single element with degree along xi,eta and neta begin linear. It contains the element indices 
            No of columns represent number of elements
            No of rows  represent number of nodes within the element

    Expected result : Front face
                        2--3
                        |  |
                        0--1  
                      Back face
                        6--7
                        |  |
                        4--5
                        control point assembly :[0,1,2,3,4,5,6,7]
    Test command : pytest test_geometry.py::test__single_element_Assembly

    Remarks : test case passed successfully
    '''
    #Input parametere are intialized
    length=1
    height=1
    width=1
    nx=2
    ny=2
    nz=2


    XI_DEGREE = 1
    ETA_DEGREE = 1
    NETA_DEGREE = 1

    N = nx
    P = ny
    Q = nz
    #Base on input parameters,the knot vector and knot span are generated
    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS = C.crtpts_coordinates()

    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()

    #Control point assembly is build based on the parameters obtained from INPUT class
    element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                            ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    
    expected_output=array([[0, 1, 2, 3, 4, 5, 6, 7]])
    number_of_elements=len(element_indicies[:,0])
    assert (number_of_elements==1 and array_equiv(expected_output,element_indicies)) is True


def test__element_Assembly_C0_continuity():
    '''
    UNIT TESTING
    Aim: Control point assembly is tested for knot vector and degree along xi,eta and neta begin linear. It contains the element indices 
            No of columns represent number of elements
            No of rows  represent number of nodes within the element

    Expected result :
    back face
    3--4--5
    |  |  |
    0--1--2

    front face
    9--10--11
    |   |   |
    6---7---8

    element_1= 0,1,3,4,6,7,9,10
    element_2=1,2,4,5,7,8,10,11

    Test command : pytest test_geometry.py::test__element_Assembly_C0_continuity

    Remarks : test case passed successfully

    '''
    #Inputs parameters are intialized
    length=1
    height=1
    width=1
    nx=3
    ny=2
    nz=2


    XI_DEGREE = 1
    ETA_DEGREE = 1
    NETA_DEGREE = 1

    N = nx
    P = ny
    Q = nz
    #Knot vector and span are obtained from the input parameters
    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS = C.crtpts_coordinates()


    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()

    #Control point assembly is build based on the parameters obtained from INPUT class
    element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                            ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    
    expected_output=array([[ 0,  1,  3,  4,  6,  7,  9, 10],
                            [ 1,  2,  4,  5,  7,  8, 10, 11]])
    number_of_elements=len(element_indicies[:,0])
    no_nodes=len(element_indicies[0,:])
    expected_nodes=(XI_DEGREE+1)*(ETA_DEGREE+1)*(NETA_DEGREE+1)
    assert (number_of_elements==2 and no_nodes==expected_nodes and array_equiv(expected_output,element_indicies)) is True


def test__single_element_Assembly_C1_continuity_along_x():
    '''
    UNIT TESTING
    Aim: Control point assembly is tested for knot vector and degree along xi begin quadratic and along eta,neta begin linear. It contains the element indices 
            No of columns represent number of elements
            No of rows  represent number of nodes within the element

    Expected result :
    BACK FACE
    3--4--5
    |     |
    0--1--2

    Front face
    9--10--11
    |       |
    6---7---8
    Control point assembly: [[ 0,  1,  2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11]]
    one element having 12 nodes
    along x direction the element has C1 continuity
    along y and z having C0 continuity

    Test command : pytest test_geometry.py::test__single_element_Assembly_C1_continuity_along_x

    Remarks : test case passed successfully
    '''
    #Inputs parameters are intialized
    length=1
    height=1
    width=1
    nx=3
    ny=2
    nz=2


    XI_DEGREE = 2
    ETA_DEGREE = 1
    NETA_DEGREE = 1

    N = nx
    P = ny
    Q = nz
    #Knot vector and span are obtained from the input parameters
    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS = C.crtpts_coordinates()


    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()

    #Control point assembly is build based on the parameters obtained from INPUT class
    element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                            ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    
    expected_output=array([[ 0,  1,  2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11]])
    number_of_elements=len(element_indicies[:,0])
    no_nodes=len(element_indicies[0,:])
    expected_nodes=(XI_DEGREE+1)*(ETA_DEGREE+1)*(NETA_DEGREE+1)
    assert (number_of_elements==1 and no_nodes==expected_nodes and array_equiv(expected_output,element_indicies)) is True
