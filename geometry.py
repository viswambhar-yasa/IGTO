import numpy as np 
import pytest 

def knot_index(degree,U,knotvector):
    '''
    This function return the knot span or the index of the region in which the value(U) lies in knotvector
    INPUTS :
        Degree      :int
                order of B-splines  0-constant, 1-linear, 2-quadratic, 3-cubic
        U           - int
                The value whose index has to be found ex-0.75
        knotvector  -  array or list
                list of knot values. Ex- [0,0,0,0.5,1,1,1]
    Returns :
        index : integer
                The position of the value U in knotvector

    '''
    n=np.size(knotvector)-1
    if U in knotvector[(degree+1):(n-degree)]:
        for i, j in enumerate(knotvector):
            if j == U:
                return i
    else:    
        if degree>0:
            n=np.size(knotvector)
            index=n-degree-1
            if U == knotvector[index+1]:
                return (index)
            if U == knotvector[degree]:
                return (degree)
            for i,j in enumerate(knotvector):
                if U>j and U<knotvector[i+1]:
                    return(i)



def bspline_basis(knot_index,degree,U,knotvector):
    '''
    modified version based on Algorthim 2.2 THE NURBS book pg70

    Parameters
    ----------
    knotindex : int 
        DESCRIPTION. The default is ''.

    degree : int
        order of B-splines  0-constant, 1-linear, 2-quadratic, 3-cubic
    U : int
            The value whose basis are to be found
    knotvector : array or list
            list of knot values.
    Returns
    -------
    Basis array of dimension degree+1
        It contains non-zero cox de boor based basis function
            Ex= if degree=2 and knotindex=4 
            basis=[N 2_0, N 3_1, N 4_2]

    '''
    if degree >0 :
        alpha=np.zeros(degree+1)
        beta=np.zeros(degree+1)
        Basis=np.zeros(degree+1)
        Basis[0]=1.0 
        for i in range(1,degree+1):
            value=0
            alpha[i]=U-knotvector[(knot_index-i+1)]
            beta[i]=knotvector[knot_index+i]-U  
            j=0
            while j<i:
                if (beta[j+1]+alpha[i-j])==0:
                     temp=0
                     Basis[j] = value+(beta[j+1]*temp)
                     value = alpha[i-j]*temp
                     Basis[i] = value
                     j=j+1
                else:
                    temp= Basis[j]/(beta[j+1]+alpha[i-j])
                    Basis[j] = value+(beta[j+1]*temp)
                    value = alpha[i-j]*temp
                    Basis[i] = value
                    j=j+1
    else:
        Basis=np.zeros(degree+1)
        if U<=knotvector[knot_index+1] and U> knotvector[knot_index]:
            Basis[0]=1
        else:
            Basis[0]=0
    return Basis



def derbspline_basis(knot_index,degree,U,knotvector):
    '''
     Modified version based on the alogorithm A2.3 in NURBS Book page no.72
    Parameters
    ----------
    knot_index : integer
                The position of the value U in knotvector
    degree : int
        order of B-splines  0-constant, 1-linear, 2-quadratic, 3-cubic
    U : int
            The value whose basis are to be found
    knotvector : array or list
            list of knot values.

    Returns
    -------
    Array 
        Derivatives of Basis array of dimension degree+1
        It contains non-zero cox de boor based basis function
            Ex= if degree=2 and knotindex=4 
            der_basis=[N' 2_0, N' 3_1, N' 4_2].

    '''
    if degree>=0:
        derBasis=np.zeros(degree+1)
        if  knot_index<degree or knot_index>(int(len(knotvector)-degree)):
            derBasis=np.zeros(degree+1)
        else:
            alpha=np.zeros(degree+1)
            beta=np.zeros(degree+1)
            Basis=np.zeros((degree+1,degree+1))
            if U>knotvector[degree] and U<=knotvector[len(knotvector)-degree]:
                Basis[0,0]=1.0
            else:
                Basis[0,0]=0
            for i in range(1,degree+1):
                alpha[i]=U-knotvector[knot_index-i+1]
                beta[i]=knotvector[knot_index+i]-U 
                saved=0.0
                for j in range(0,i):
                    Basis[i,j]=beta[j+1]+alpha[i-j]    
                    temp= Basis[j,i-1]/(beta[j+1]+alpha[i-j])
                    Basis[j,i] = saved+(beta[j+1]*temp)
                    saved = alpha[i-j]*temp
                Basis[i,i] = saved  
        for r in range(0,degree+1):
            if (r>=1):
                derBasis[r]=Basis[r-1,degree-1]/Basis[degree,r-1]
            if (r<=degree-1):
                derBasis[r]=derBasis[r]-(Basis[r,degree-1]/Basis[degree,r])
        derBasis *=degree
    return np.array(derBasis)

