import numpy as np
from geometry import controlpointassembly,knot_connectivity,trilinear_der
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
    if p>=5:
        p=4
    else:
        p=p+1
    if q>=5:
        q=q+1
    else:
        q=q+1
    if r>=5:
        r=4
    else:
        r=r+1
    
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



def jacobian(Xi,Eta,Nta,Uspan,Vspan,Wspan,X,Y,Z,weights,xdegree,xknot_vector,ydegree,yknot_vector,zdegree,zknot_vector):
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
    pts=np.array([X,Y,Z])
    [dR_dxi,dR_deta,dR_dnta,nurbs]=trilinear_der(Xi,Eta,Nta,weights,xdegree,xknot_vector,ydegree,yknot_vector,zdegree,zknot_vector)
    DR_DXi=np.array([dR_dxi,dR_deta,dR_dnta])
    #jacobian12=pts@ np.transpose(DR_DXi)
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
    #jacobian1=jacobian12
    det_jacobian1=np.linalg.det(jacobian1)

    #print(jacobian1)
    return jacobian1,det_jacobian1,det_jacobian2,DR_DXi,nurbs



def strain_displacement(Xi,Eta,Neta,jacobian1,det_jacobian1,DR_DXi,nurbs,nn):
    '''
    

    Parameters
    ----------
    Xi : TYPE
        DESCRIPTION.
    Eta : TYPE
        DESCRIPTION.
    Neta : TYPE
        DESCRIPTION.
    jacobian1 : TYPE
        DESCRIPTION.
    det_jacobian1 : TYPE
        DESCRIPTION.
    DR_DXi : TYPE
        DESCRIPTION.
    NURBS : TYPE
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    B : TYPE
        DESCRIPTION.
    R : TYPE
        DESCRIPTION.
    '''
    
    #if det_jacobian1==0:
    #    raise Exception("Jacobian cannot be zero")
    inverse_jacobian1=np.linalg.inv(jacobian1)
    DR_DX=np.matmul(inverse_jacobian1,DR_DXi)

    B=np.zeros((6,int(3*np.size(DR_DX,1))))
    R=np.zeros((3,int(3*np.size(DR_DX,1))))
    for i in range(int(np.size(DR_DX,1))):
        x=int(3*i)
        y=int(3*i+1)
        z=int(3*i+2)
        B[0,x]=DR_DX[0,i]
        B[1,y]=DR_DX[1,i]
        B[2,z]=DR_DX[2,i]
        B[3,x]=DR_DX[1,i]
        B[3,y]=DR_DX[0,i]
        B[4,y]=DR_DX[2,i]
        B[4,z]=DR_DX[1,i]
        B[5,x]=DR_DX[2,i]
        B[5,z]=DR_DX[0,i]
        R[0,x]=nurbs[i]
        R[0,y]=nurbs[i]
        R[2,z]=nurbs[i]
    #print(B)
    return B,R


def Compliance_matrix(E=100000,v=0.3):
    '''
    

    Parameters
    ----------
    E : TYPE
        DESCRIPTION.
    v : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    C=np.zeros((6,6))
    C[0,0]=C[1,1]=C[2,2]=(E/((1+v)*(1-2*v)))*(1-v)
    C[3,3]=C[4,4]=C[5,5]=E/((1+v)*1)
    C[0,1]=C[0,2]=C[1,0]=C[2,0]=C[1,2]=C[2,1]=(E/((1+v)*(1-2*v)))*v
    return C

def element_routine(X,Y,Z,weights,E,v,Uspan,Vspan,Wspan,xdegree,xknot_vector,ydegree,yknot_vector,zdegree,zknot_vector):
    Gauss_points,gauss_weights=gauss_quadrature(xdegree,ydegree,zdegree)
    C=Compliance_matrix(E,v)
    Ke=np.zeros((3*((xdegree+1)*(ydegree+1)*(zdegree+1)),3*((xdegree+1)*(ydegree+1)*(zdegree+1))))
    nn=((xdegree+1)*(ydegree+1)*(zdegree+1))
    for i in range(len(gauss_weights)):
        gp=Gauss_points[i,:]
        Wt=gauss_weights[i]
        Xi=unittoparametric(gp[0],Uspan)
        Eta=unittoparametric(gp[1],Vspan)
        Neta=unittoparametric(gp[2],Wspan)
        [jacobian1,det_jacobian1,det_jacobian2,DR_DXi,nurbs]=jacobian(Xi,Eta,Neta,Uspan,Vspan,Wspan,X,Y,Z,weights,xdegree,xknot_vector,ydegree,yknot_vector,zdegree,zknot_vector)
        B,R=strain_displacement(Xi,Eta,Neta,jacobian1,det_jacobian1,DR_DXi,nurbs,nn)
        temp=np.transpose(B)@C
        Ke+=temp@B*(det_jacobian1*det_jacobian2*Wt)
    return Ke,nurbs,R


def assemble(K_G,K_E,elindices,ncp,K_disp=False):
    '''
    

    Parameters
    ----------
    K_G : TYPE
        DESCRIPTION.
    K_E : TYPE
        DESCRIPTION.
    elindices : TYPE
        DESCRIPTION.
    ncp : TYPE
        DESCRIPTION.

    Returns
    -------
    K_G : TYPE
        DESCRIPTION.

    '''
    
    elindices=np.array(elindices)
    dof=3
    gindices=np.sort(np.concatenate((elindices*dof,dof*elindices+1,dof*elindices+2)))
    for i,ii in enumerate(gindices):
        for j,jj in enumerate(gindices):
            K_G[ii,jj]=K_G[ii,jj]+K_E[i,j]
    if K_disp:        
        print('Building Global Stiffness matrix \n')
    return K_G

def apply_BC(K_G,F_E,fixed_dof,load_dof,P,abc_disp=False):
    '''
    

    Parameters
    ----------
    K_G : TYPE
        DESCRIPTION.
    F_E : TYPE
        DESCRIPTION.
    fixed_dof : TYPE
        DESCRIPTION.
    load_dof : TYPE
        DESCRIPTION.
    P : TYPE
        DESCRIPTION.

    Returns
    -------
    reduced_GK : TYPE
        DESCRIPTION.
    reduced_F_E : TYPE
        DESCRIPTION.

    '''
    #reduced_GK=np.delete(K_G,fixed_dof,axis=0) #axis=0 is row
    #reduced_GK=np.delete(reduced_GK,fixed_dof,axis=1) 
    reduced_GK=np.delete(np.delete(K_G, fixed_dof, 0),fixed_dof , 1)
    #axis=0 is row
    F_E[load_dof]=P
    reduced_F_E=np.delete(F_E,fixed_dof,axis=0)
    if abc_disp:
        print('Applying Boundary Conditions \n')
    return reduced_GK,reduced_F_E

#print(gauss_quadrature(1,1,1))