import numpy as np
from geometry import trilinear_der
from Inputs import *
import pytest

def gauss_quadrature(p,q,r):
    '''
    

    Parameters
    ----------
    p : TYPE
        DESCRIPTION.
    q : TYPE
        DESCRIPTION.
    r : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    '''
    [xgauss_points,xweights]=np.polynomial.legendre.leggauss(p)
    [ygauss_points,yweights]=np.polynomial.legendre.leggauss(q)
    [zgauss_points,zweights]=np.polynomial.legendre.leggauss(r)
    gausspoints=np.zeros((p*q*r,3))
    weights=np.zeros(p*q*r)
    n=0
    for i in range(r):
        for j in range(q):
            for k in range(p):
                gausspoints[n,:]=[xgauss_points[k],ygauss_points[j],zgauss_points[i]]
                weights[n]=xweights[k]*yweights[j]*zweights[i]
                n+=1
    return gausspoints,weights


