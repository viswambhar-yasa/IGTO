from geometry import knot_index,bspline_basis,derbspline_basis,trilinear_der
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