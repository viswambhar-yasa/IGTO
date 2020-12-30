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
            return derBasis
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


def trilinear_der(Ux,Uy,Uz,weights,xdegree,xknotvector,ydegree,yknotvector,zdegree,zknotvector):
    '''
    

    Parameters
    ----------
    Ux : TYPE
        DESCRIPTION.
    Uy : TYPE
        DESCRIPTION.
    Uz : TYPE
        DESCRIPTION.
    weights : TYPE
        DESCRIPTION.
    xdegree : TYPE, optional
        DESCRIPTION. The default is XI_DEGREE.
    xknotvector : TYPE, optional
        DESCRIPTION. The default is XI_KNOTVECTOR.
    ydegree : TYPE, optional
        DESCRIPTION. The default is ETA_DEGREE.
    yknotvector : TYPE, optional
        DESCRIPTION. The default is ETA_KNOTVECTOR.
    zdegree : TYPE, optional
        DESCRIPTION. The default is ETA_DEGREE.
    zknotvector : TYPE, optional
        DESCRIPTION. The default is ETA_KNOTVECTOR.

    Returns
    -------
    None.

    '''
    x_index=knot_index(xdegree,Ux,xknotvector)
    y_index=knot_index(ydegree,Uy,yknotvector)
    z_index=knot_index(zdegree,Uz,zknotvector)


    Nx=bspline_basis(x_index,xdegree,Ux,xknotvector)
    Ny=bspline_basis(y_index,ydegree,Uy,yknotvector)
    Nz=bspline_basis(z_index,zdegree,Uz,zknotvector)


    DNx=derbspline_basis(x_index,xdegree,Ux,xknotvector)
    DNy=derbspline_basis(y_index,ydegree,Uy,yknotvector)
    DNz=derbspline_basis(z_index,zdegree,Uz,zknotvector)


    W=0
    W_dx=0
    W_dy=0
    W_dz=0
    p=0
    windex=0

    for k in range(zdegree+1):
        for j in range(ydegree+1):
            for i in range(xdegree+1):
                W+=Nx[i]*Ny[j]*Nz[k]*weights[windex]
                W_dx+=DNx[i]*Ny[j]*Nz[k]*weights[windex]
                W_dy+=Nx[i]*DNy[j]*Nz[k]*weights[windex]
                W_dz+=Nx[i]*Ny[j]*DNz[k]*weights[windex]
                windex+=1

    dR_dx=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
    dR_dy=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
    dR_dz=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
    R=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))

    w=0
    windex=0

    if W==0:
        dR_dx=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
        dR_dy=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
        dR_dz=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
        R=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
    else:
        for k in range(zdegree+1):
            for j in range(ydegree+1):
                for i in range(xdegree+1):
                    w=weights[windex]/(W*W)
                    R[p]=(Nx[i]*Ny[j]*Nz[k]*W*w)
                    dR_dx[p]=(DNx[i]*Ny[j]*Nz[k]*W-Nx[i]*Ny[j]*Nz[k]*W_dx)*w
                    dR_dy[p]=(Nx[i]*DNy[j]*Nz[k]*W-Nx[i]*Ny[j]*Nz[k]*W_dy)*w
                    dR_dz[p]=(Nx[i]*Ny[j]*DNz[k]*W-Nx[i]*Ny[j]*Nz[k]*W_dz)*w
                    p=p+1
                    windex+=1

    return [dR_dx,dR_dy,dR_dz,R]


def elementorder(numx,numy,numz):
    '''
    

    Parameters
    ----------
    numx : TYPE
        DESCRIPTION.
    numy : TYPE
        DESCRIPTION.
    numz : TYPE
        DESCRIPTION.

    Returns
    -------
    el_order : TYPE
        DESCRIPTION.

    '''
    index=0
    el_order=np.zeros((numx,numz,numy))
    for i in range(numz):
        for j in range(numy):
            for k in range(numx):
                el_order[k,i,j]=index
                index+=1
    return el_order

def elementorder1(numx,numy,numz):
    '''
    

    Parameters
    ----------
    numx : TYPE
        DESCRIPTION.
    numy : TYPE
        DESCRIPTION.
    numz : TYPE
        DESCRIPTION.

    Returns
    -------
    el_order : TYPE
        DESCRIPTION.

    '''
    index=0
    el_order=np.zeros((numx,numz,numy))
    for k in range(numz):
        for j in range(numy):
            for i in range(numx):
                el_order[i,k,j]=index
                index+=1
    return el_order

#print(elementorder1(2,2,2))


def knot_connectivity(n,p,q,knotconnectivityU,knotconnectivityV,knotconnectivityW):
    '''
    

    Parameters
    ----------
    element : TYPE
        DESCRIPTION.
    n : TYPE, optional
        DESCRIPTION. The default is N.
    p : TYPE, optional
        DESCRIPTION. The default is P.
    q : TYPE, optional
        DESCRIPTION. The default is Q.
    knotconnectivityU : TYPE, optional
        DESCRIPTION. The default is XI_KNOTCONNECTIVITY.
    knotconnectivityV : TYPE, optional
        DESCRIPTION. The default is ETA_KNOTCONNECTIVITY.
    knotconnectivityW : TYPE, optional
        DESCRIPTION. The default is NETA_KNOTCONNECTIVITY.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    index=np.zeros((n*p*q,3))
    count=0
    for k in range(len(knotconnectivityW)):
        for j in range(len(knotconnectivityV)):
            for i in range(len(knotconnectivityU)):
                index[count,:]=[i,j,k]
                count+=1
    index=index.astype(int)
    return  index


def controlpointassembly(n,p,q,nU,nV,nW,xdegree,ydegree,zdegree,knotconnectivityU,knotconnectivityV,knotconnectivityW):#modigy this 
    '''
    for 2 elements
    
    front face numbering
    3-4-5
    0-1-2

    back face numbering
    9-10-11
    6-7-8
    '''
    elements_assembly=np.zeros(((nU*nV*nW),(xdegree+1)*(ydegree+1)*(zdegree+1)))
    elements_order=elementorder(n,p,q)
    a=0 
    for i in range(nW):
        for j in range(nV):
            for k in range(nU):
                c=0
                for l in range(len(knotconnectivityW[i,:])):
                    for m in range(len(knotconnectivityV[j,:])):
                        for n in range(len(knotconnectivityU[k,:])):
                            elements_assembly[a,c]=elements_order[knotconnectivityU[k,n],knotconnectivityW[i,l], knotconnectivityV[j,m]]
                            c+=1
                a+=1         
    elements_assembly=elements_assembly.astype(int)
    return elements_assembly

'''
knotvector=[0, 0, 0, 0,0, 1/4, 1/3,1/2, 3/4,1, 1, 1, 1, 1 ]
degree=4
U=0.6
knotindex=knot_index(degree,U,knotvector)
print(derbspline_basis(knotindex,degree,U,knotvector))
'''