from element_routine import gauss_quadrature
from numpy import array,array_equiv
import pytest

def test__gauss_quadrature_1point_true():
    gausspoints,weights=gauss_quadrature(1,1,1)
    assert( array_equiv(gausspoints,array([0.0,0.0,0.0])) and array_equiv(weights,array([8.0])) )is True
