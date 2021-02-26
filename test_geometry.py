from Preprocessing import Inputs
from geometry import knot_index,bspline_basis,derbspline_basis,trilinear_der,controlpointassembly
from numpy import array,array_equiv,sum,equal,round,ones,float,zeros,all
import pytest

def test__knot_index_true():
    '''
    Values obtained for NURBS book
    '''
    knotvector=[0,0,0,0.14285714,0.28571429,0.42857143,0.57142857,0.71428571,0.85714286,1,1,1]
    degree=2
    U=0.3
    assert (knot_index(degree,U,knotvector)==4) is True

def test__Bspline_basis_equal_true():
    '''
    Values obtained for NURBS book
    '''
    knotvector=[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
    degree=2
    U=5/2
    knotindex=knot_index(degree,U,knotvector)
    output=bspline_basis(knotindex,degree,U,knotvector)
    assert (array_equiv(output,[0.125,0.75,0.125])) is True

def test__Bspline_basis_sum_equal_to_one_true():
    '''
    Values obtained for NURBS book
    '''
    knotvector=[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
    degree=2
    U=5/3
    knotindex=knot_index(degree,U,knotvector)
    output=sum(bspline_basis(knotindex,degree,U,knotvector))
    assert (int(output)==1) is True

def test__Bspline_basis_all_values_positive_true():
    '''
    Values obtained for NURBS book
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



def test__derbspline_basis_sum_equal_zero_true():
    '''#
    Values obtained for NURBS book
    '''
    knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
    degree=3
    U=0.6
    knotindex=knot_index(degree,U,knotvector)
    output=round(sum(derbspline_basis(knotindex,degree,U,knotvector)),8)
    assert (float(output)==0.0) is True


def test__derbspline_basis_equal_true():
    '''
    Values obtained for NURBS book
    '''
    knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
    degree=3
    U=0.6
    knotindex=knot_index(degree,U,knotvector)
    output=round(derbspline_basis(knotindex,degree,U,knotvector),2 )
    expected_output=array([-0.72, -2.24,  2.48,  0.48])
    assert (array_equiv(output,expected_output)) is True

def test__trilinear_der_Basis_sum_equal_to_one_true():
    '''
    Values obtained from NURBS book
    '''
    WEIGHTS=ones(8)
    Xi_degree=Eta_degree=Neta_degree=1
    Xi_knotvector=Eta_knotvector=Eta_knotvector=[0,0,0.5,1,1]
    DR_DX,DR_DY,DR_DZ,R =trilinear_der(0.6,0.5,0.3,WEIGHTS,Xi_degree,Xi_knotvector,Eta_degree,Eta_knotvector,Neta_degree,Eta_knotvector)
    output=sum(R)
    assert (float(output)==1.0) is True


def test__trilinear_der_Basis_less_than_zero_false():
    
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
    Values obtained from NURBS book
    '''
    Xi_degree=Eta_degree=Neta_degree=2
    WEIGHTS=ones((Xi_degree+1)**3)
    Xi_knotvector=Eta_knotvector=Eta_knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
    DR_DXI,DR_DETA,DR_DNETA,R =trilinear_der(1/4,2/3,3/5,WEIGHTS,Xi_degree,Xi_knotvector,Eta_degree,Eta_knotvector,Neta_degree,Eta_knotvector)
    output=float(round(sum(DR_DXI),10))
    assert (output==0.0) is True


def test__trilinear_der_ETA_sum_equal_to_zero_true():
    '''
    Values obtained from NURBS book
    '''
    Xi_degree=Eta_degree=Neta_degree=2
    WEIGHTS=ones((Xi_degree+1)**3)
    Xi_knotvector=Eta_knotvector=Eta_knotvector=[0,0,0,0.16666667,0.33333333, 0.5,0.66666667, 0.83333333,1,1,1]
    DR_DXI,DR_DETA,DR_DNETA,R =trilinear_der(0.5,0.8,0.2,WEIGHTS,Xi_degree,Xi_knotvector,Eta_degree,Eta_knotvector,Neta_degree,Eta_knotvector)
    output=float(round(sum(DR_DETA),10))
    assert(output==0) is True

def test__trilinear_der_NETA_sum_equal_to_zero_true():
    '''
    Values obtained from NURBS book
    '''
    Xi_degree=Eta_degree=Neta_degree=3
    WEIGHTS=ones((Xi_degree+1)**3)
    Xi_knotvector=Eta_knotvector=Eta_knotvector=[0,0,0,0,0.33333333, 0.66666667,1,1,1,1]
    DR_DXI,DR_DETA,DR_DNETA,R =trilinear_der(0.5,0.8,0.2,WEIGHTS,Xi_degree,Xi_knotvector,Eta_degree,Eta_knotvector,Neta_degree,Eta_knotvector)
    output=float(round(sum(DR_DNETA),10))
    assert(output==0) is True


def test__single_element_Assembly():
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

    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS = C.crtpts_coordinates()

    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()


    element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                            ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    
    expected_output=array([[0, 1, 2, 3, 4, 5, 6, 7]])
    number_of_elements=len(element_indicies[:,0])
    assert (number_of_elements==1 and array_equiv(expected_output,element_indicies)) is True


def test__element_Assembly_C0_continuity():
    '''
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
    '''
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

    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS = C.crtpts_coordinates()


    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()


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
    BACK FACE
    3--4--5
    |     |
    0--1--2

    Front face
    9--10--11
    |       |
    6---7---8

    one element having 12 nodes
    along x direction the element has C1 continuity
    along y and z having C0 continuity
    '''
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

    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS = C.crtpts_coordinates()


    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()


    element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                            ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    
    expected_output=array([[ 0,  1,  2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11]])
    number_of_elements=len(element_indicies[:,0])
    no_nodes=len(element_indicies[0,:])
    expected_nodes=(XI_DEGREE+1)*(ETA_DEGREE+1)*(NETA_DEGREE+1)
    assert (number_of_elements==1 and no_nodes==expected_nodes and array_equiv(expected_output,element_indicies)) is True
