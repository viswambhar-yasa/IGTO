import numpy as np
from geometry import controlpointassembly,knot_connectivity,trilinear_der
from Inputs import *
import pytest

def gauss_quadrature(p=XI_DEGREE,q=ETA_DEGREE,r=NETA_DEGREE):
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
    if p>=5:
        p=4
    else:
        p=p+1
    if q>=5:
        q=4
    else:
        q=q+1
    if r>=5:
        r=4
    else:
        r=p+1
    
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



def unittoparametric(gauss_point,span):
    '''
    

    Parameters
    ----------
    gauss_point : TYPE
        DESCRIPTION.
    span : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return (np.dot(0.5,((span[1]-span[0])*gauss_point+(span[1]+span[0]))))



def jacobian(Xi,Eta,Nta,Uspan,Vspan,Wspan,X,Y,Z,weights,xdegree=XI_DEGREE,xknot_vector=XI_KNOTVECTOR,ydegree=ETA_DEGREE,yknot_vector=ETA_KNOTVECTOR,zdegree=NETA_DEGREE,zknot_vector=NETA_KNOTVECTOR):
    '''
    

    Parameters
    ----------
    Xi : TYPE
        DESCRIPTION.
    Eta : TYPE
        DESCRIPTION.
    Nta : TYPE
        DESCRIPTION.
    Uspan : TYPE
        DESCRIPTION.
    Vspan : TYPE
        DESCRIPTION.
    Wspan : TYPE
        DESCRIPTION.
    X : TYPE
        DESCRIPTION.
    Y : TYPE
        DESCRIPTION.
    Z : TYPE
        DESCRIPTION.
     : TYPE
        DESCRIPTION.
    weights : TYPE
        DESCRIPTION.
    xdegree : TYPE, optional
        DESCRIPTION. The default is XI_DEGREE.
    xknot_vector : TYPE, optional
        DESCRIPTION. The default is XI_KNOTVECTOR.
    ydegree : TYPE, optional
        DESCRIPTION. The default is ETA_DEGREE.
    yknot_vector : TYPE, optional
        DESCRIPTION. The default is ETA_KNOTVECTOR.
    zdegree : TYPE, optional
        DESCRIPTION. The default is NETA_DEGREE.
    zknot_vector : TYPE, optional
        DESCRIPTION. The default is NETA_KNOTVECTOR.

    Returns
    -------
    list
        DESCRIPTION.

    '''
    DXi_dxi=0.5*(Uspan[1]-Uspan[0])
    DEta_dEta=0.5*(Vspan[1]-Vspan[0])
    DNta_dNta=0.5*(Wspan[1]-Wspan[0])
    det_jacobian2=DXi_dxi*DEta_dEta*DNta_dNta

    [dR_dxi,dR_deta,dR_dnta,nurbs]=trilinear_der(Xi,Eta,Nta,weights)
    DR_DXi=np.array([dR_dxi,dR_deta,dR_dnta])

    jacobian1=np.zeros((3,3))

    dRxi_dx=np.dot(dR_dxi,X)
    jacobian1[0,0]=dRxi_dx
    
    dRxi_dy=np.dot(dR_dxi,Y)
    jacobian1[0,1]=dRxi_dy
    
    dRxi_dz=np.dot(dR_dxi,Z)
    jacobian1[0,2]=dRxi_dz
    
    dReta_dx=np.dot(dR_deta,X)
    jacobian1[1,0]=dReta_dx
    
    dReta_dy=np.dot(dR_deta,Y)
    jacobian1[1,1]=dReta_dy
    
    dReta_dz=np.dot(dR_deta,Z)
    jacobian1[1,2]=dReta_dz
    
    dRNta_dx=np.dot(dR_dnta,X)
    jacobian1[2,0]=dRNta_dx
    
    dRNta_dy=np.dot(dR_dnta,Y)
    jacobian1[2,1]=dRNta_dy
    
    dRNta_dz=np.dot(dR_dnta,Z)
    jacobian1[2,2]=dRNta_dz

    det_jacobian1=np.linalg.det(jacobian1)


    return [jacobian1,det_jacobian1,det_jacobian2,DR_DXi,nurbs]


