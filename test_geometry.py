from geometry import knot_index,bspline_basis,derbspline_basis
from numpy import array, array_equiv,sum,not_equal
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
    output=sum(derbspline_basis(knotindex,degree,U,knotvector))
    assert (equal( output, 0)) is True


def test__derbspline_basis_equal_true():
    '''
    Values obtained for NURBS book
    '''
    knotvector=[0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1 ]
    degree=3
    U=0.6
    knotindex=knot_index(degree,U,knotvector)
    output=sum(derbspline_basis(knotindex,degree,U,knotvector))
    assert (array_equiv(output,array([-0.72, -2.24,  2.48,  0.48]))) is False

