from element_routine import gauss_quadrature,unittoparametric,jacobian,assemble,element_routine,apply_BC 
from numpy import array,array_equiv,zeros,linalg,sin,cos,radians,transpose
from Inputs import *
from geometry import knot_connectivity,controlpointassembly
from boundary_conditions import BC_Switcher
import pytest

def test__gauss_quadrature_1point_true():
    gausspoints,weights=gauss_quadrature(0,0,0)
    assert( array_equiv(gausspoints,array([0.0,0.0,0.0])) and array_equiv(weights,array([8.0])) )is True

def test__unit_to_parametric_space_true():
    '''
    value calculated anaytically 

    '''
    output=float(round(unittoparametric(-0.5773,[0,0.2]),5))
    assert (output==float(0.04227)) is True



def test_Jacobian_patch_testing_rotation():
    '''
    rotation of the jacbian matrix using Bunge convenstion should change the determinat as it is volume conserving
    '''
    
    X=unittoparametric(-0.5773,[0,1])
    E=unittoparametric(0.5773,[0,1])
    N=unittoparametric(-0.5773,[0,1])
    Wspa=Vspa=Uspa=array([0,1])
    control_points=array([[0., 0., 0. ,1.],
                    [1., 0., 0. ,1.],
                    [0., 1., 0., 1.],
                    [1., 1., 0., 1.],
                    [0., 0., 1., 1.],
                    [1., 0., 1., 1.],
                    [0., 1., 1., 1.],
                    [1., 1., 1., 1.]])
    xx=control_points[:,0]
    yy=control_points[:,1]
    zz=control_points[:,2]
    weigt=control_points[:,3]
    xdef=1
    XKV=[0, 0, 1, 1]
    jacobian1,N1,N2,N3,N4=jacobian(X,E,N,Uspa,Vspa,Wspa,xx,yy,zz,weigt,xdef,XKV,xdef,XKV,xdef,XKV)
    output=zeros(360)
    for i in range(0,360):
        a=radians(i)
        Q1=array([[cos(a),sin(a) ,0],[-sin(a),cos(a) ,0],[0,0,1]])
        Q2=array([[cos(a),0,-sin(a)],[0,1,0],[sin(a),0,cos(a)]])
        Q3=array([[cos(a),sin(a) ,0],[-sin(a),cos(a) ,0],[0,0,1]])
        Q=Q1@Q2@Q3
        Q_T=transpose(Q)
        output[i]=round(linalg.det(Q_T@jacobian1@Q),10)
    actual_det_jacobian=round(linalg.det(jacobian1),10)
    assert all(output==actual_det_jacobian) is True  

