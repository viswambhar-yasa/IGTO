#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#-------------------------------------------------------------------------------------------#
#GEOMETRY - Python file used to generate geometry parameter require to build element routine
#-------------------------------------------------------------------------------------------#


import numpy as np 
import pytest 


def knot_index(degree,U,knotvector):
    '''
    This function return the knot span or the index of the region in which the value(U) lies in knotvector.

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
    
    Test Case:
    Test command - pytest test_geometry.py::test__knot_index_true

    '''
    
    n=np.size(knotvector)-1   #length of the knotvector

    # knot index is calculated using simple search method
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
    Implementation of equation 2.2 for the documentation

    Generates Basis functions.

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
    Basis: Array dimension degree+1   
            It contains non-zero cox de boor based basis function
                Ex= if degree=2 and knotindex=4 
                    basis=[N 2_0, N 3_1, N 4_2]

    Test case
    ---------
    test commands - 
                    pytest test_geometry.py::test__Bspline_basis_equal_true
                    
                    pytest test_geometry.py::test__Bspline_basis_sum_equal_to_one_true

                    pytest test_geometry.py::test__Bspline_basis_all_values_positive_true
                    
    '''
    if degree >0 :
        #initization of the base vectors 
        alpha=np.zeros(degree+1)
        beta=np.zeros(degree+1)
        Basis=np.zeros(degree+1)
        Basis[0]=1.0 
       #higher order Basis are calculated as in figure 5
        for i in range(1,degree+1):
            value=0
            #the range of the basis is calculated and stored in alpha and beta
            alpha[i]=U-knotvector[(knot_index-i+1)]
            beta[i]=knotvector[knot_index+i]-U  
            j=0
            #implementation of equ (2.2)  
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
    #This section computes lower order basis function
    else:
        Basis=np.zeros(degree+1)
        if U<=knotvector[knot_index+1] and U> knotvector[knot_index]:
            Basis[0]=1
        else:
            Basis[0]=0
    return Basis


def derbspline_basis(i,p,xi,XI,g=1):#modify
    """
    Modified version based on the alogorithm A2.3 in NURBS Book page no.72
    Generates derivatives of basis function.

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
      """
    #Inititation of dimentions for 2D matrices
    ndu=np.zeros((p+1,p+1))
    ders=np.zeros((p+1,p+1))
    a=np.zeros((2,p+1))
    left =np.zeros(p+2)
    right =np.zeros(p+2)
    
    ndu[0][0]=1.0
    for j in range(1,p+1):
        left[j] = xi - XI[i+1-j]
        right[j] = XI[i+j] - xi
        saved=0.0
        for r in range(j):
            #Lower triangle
            ndu[j][r] = right[r+1]+left[j-r]
            temp=ndu[r][j-1]/ndu[j][r]
            #Upper triangle
            ndu[r][j] = saved+(right[r+1]*temp)
            saved=left[j-r]*temp
        ndu[j][j] = saved
    for j in range (p+1): #Load the basis functions
        ders[0][j]=ndu[j][p]
    #This secion computes the derivatives
    for r in range(p+1):
        s1=0
        s2=1 #Alternative rows in array a
        a[0][0] = 1.0
        #Loop to compute kth derivative
        for k in range(1,g+1):
            d=0.0
            rk=r-k
            pk=p-k
            if(r>=k):
                a[s2][0]=a[s1][0]/ndu[pk+1][rk]
                d=a[s2][0]*ndu[rk][pk]
            if(rk>=-1):
                j1=1
            else:
                j1=-rk
            if(r-1<=pk):
                j2=k-1
            else:
                j2=p-r
            for j in range (j1,j2+1):
                a[s2][j] =(a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j]
                d += (a[s2][j]*ndu[rk+j][pk])
            if(r<=pk):
                a[s2][k]=-a[s1][k-1]/ndu[pk+1][r]
                d+=(a[s2][k]*ndu[r][pk])
            ders[k][r]=d
            #Switch rows
            j=s1
            s1=s2
            s2=j
            #Multiply through by the correct factors
    r=p
    for k in range(1,g+1):
        for j in range(p+1):
            ders[k][j] =ders[k][j]* r
        r =r* (p-k)
    derivatives=ders[1,:]
    return np.array(derivatives)

def trilinear_der(Ux,Uy,Uz,weights,xdegree,xknotvector,ydegree,yknotvector,zdegree,zknotvector):
    '''
    NURBS volume which is a tensor product of NURBS basis along xi,eta,neta direction.

    Parameters
    ----------

    Ux,Uy,Uz : float
                knot index along xi,eta,neta direction

    weights : Array 
                An array of NURBS weights defined for each knot.

    xdegree,ydegree,zdegree : int
                                degree of the knotvector along xi,eta,neta
                                     0-constant, 1-linear, 2-quadratic, 3-cubic

    xknotvector,yknotvector,zknotvector : Array
                                            knotvectors defined along xi,eta,neta

    Returns
    -------
    R - Array.
        Trivariant function 

    dR_dx,dR_dy,dR_dz -Array
        Derivative of Trivariant function along xi,eta,neta direction.

    Test cases
    ----------
    test command -
                    pytest test_geometry.py::test__trilinear_der_Basis_sum_equal_to_one_true

                    pytest test_geometry.py::test__trilinear_der_Basis_less_than_zero_false

                    pytest test_geometry.py::test__trilinear_der_XI_sum_equal_to_zero_true

                    pytest test_geometry.py::test__trilinear_der_ETA_sum_equal_to_zero_true

                    pytest test_geometry.py::test__trilinear_der_NETA_sum_equal_to_zero_true      
    
    '''
    #Knot index are calculated for Ux,Uy,Uz along xi,eta,neta respectively
    x_index=knot_index(xdegree,Ux,xknotvector)
    y_index=knot_index(ydegree,Uy,yknotvector)
    z_index=knot_index(zdegree,Uz,zknotvector)

    # Basis are calculated along xi,eta,neta
    Nx=bspline_basis(x_index,xdegree,Ux,xknotvector)
    Ny=bspline_basis(y_index,ydegree,Uy,yknotvector)
    Nz=bspline_basis(z_index,zdegree,Uz,zknotvector)

    #Derivtive of the Basis are calculated along xi,eta ,neta
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
                W+=Nx[i]*Ny[j]*Nz[k]*weights[windex]      #Based on equ. 3.9
                W_dx+=DNx[i]*Ny[j]*Nz[k]*weights[windex]  #Based on equ. 3.10
                W_dy+=Nx[i]*DNy[j]*Nz[k]*weights[windex]  #Based on equ. 3.11
                W_dz+=Nx[i]*Ny[j]*DNz[k]*weights[windex]  #Based on equ. 3.12
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
                    R[p]=(Nx[i]*Ny[j]*Nz[k]*W*w)                               #Based on equ. 3.5
                    dR_dx[p]=(DNx[i]*Ny[j]*Nz[k]*W-Nx[i]*Ny[j]*Nz[k]*W_dx)*w   #Based on equ. 3.6
                    dR_dy[p]=(Nx[i]*DNy[j]*Nz[k]*W-Nx[i]*Ny[j]*Nz[k]*W_dy)*w   #Based on equ. 3.7
                    dR_dz[p]=(Nx[i]*Ny[j]*DNz[k]*W-Nx[i]*Ny[j]*Nz[k]*W_dz)*w   #Based on equ. 3.8
                    p=p+1
                    windex+=1

    return [dR_dx,dR_dy,dR_dz,R]


def elementorder(numx,numy,numz):
    '''
    This function generates a 3D matrix with element ordered .This is used in the generation of control point asssembly. 

    Parameters
    ----------
    numx,numy,numz : int
                    Number of element in x,y and z direction.

    Returns
    -------
    el_order : 3D Array
                The elements are numbered and order in 3D matrix.
                Ex:
                nx,ny,nz=2,2,1
                el_order-[[[0 2]],
                            [[1 3]]]
    '''
    index=0
    #Inititation of dimensions for 3D matrices
    el_order=np.zeros((numx,numz,numy))
    for i in range(numz):
        for j in range(numy):
            for k in range(numx):
                el_order[k,i,j]=index  #element ordering in 3d matrix is performed 
                index+=1
    return el_order



def knot_connectivity(nx,ny,nz,knotconnectivityU,knotconnectivityV,knotconnectivityW):
    '''
    Based on Isogeometric analysis: an overview and computer implementation aspects paper page - 50,51 
    IGFEM 6.1. Data structure- Listing 2: - Which is implemented for a 2D element extended it for 3D elements.

    This function is used in element routine to get respective knots present within each elements.

    Parameters
    ----------

    nx,ny,nz : int
             Number of element in x,y and z direction.
    knotconnectivityU,knotconnectivityV,knotconnectivityW : Array
                                                                An array which contains the length of each element and number of knot in it.
                                                                Ex-knotvector:[0,1,2,3,4,5,6]         
    Returns
    -------
    index : 2d array
            Based on the number of element and knot span of each element in x,y and z direction .The knots spans are generated.
    
    '''
    #Inititation of dimensions for 2D matrices
    index=np.zeros((nx*ny*nz,3))
    count=0
    for k in range(len(knotconnectivityW)):
        for j in range(len(knotconnectivityV)):
            for i in range(len(knotconnectivityU)):
                index[count,:]=[i,j,k]     #knot index for each element is generated.
                count+=1
    index=index.astype(int)
    return  index


def controlpointassembly(nx,ny,nz,nU,nV,nW,xdegree,ydegree,zdegree,knotconnectivityU,knotconnectivityV,knotconnectivityW):
    '''
    Based on Isogeometric analysis: an overview and computer implementation aspects paper page - 50,51 
    IGFEM 6.1. Data structure- Listing 2: - Which is implemented for a 2D element extended it for 3D elements.

    This function generates a control point assembly matrix used in element routine 

    EX:(nx,ny,nz):(2,1,1) degree:(2,1,1)
        front face numbering
        3--4--5
        |     |
        0--1--2

        back face numbering
        9--10--11
        |      |
        6--7--8

    
    Parameters
    ----------

    nx,ny,nz : int
             Number of element in x,y and z direction.
    nU,nV,nW : int
                length of knot vector along xi,eta,neta direction.
    xdegree,ydegree,zdegree : int
                                degree of the knotvector along xi,eta,
                                 0-constant, 1-linear, 2-quadratic, 3-cubic

    knotconnectivityU,knotconnectivityV,knotconnectivityW : Array
                                                                An array which contains the length of each element and number of knot in it.
                                                                Ex-knotvector:[0,1,2,3,4,5,6]         
    Returns
    -------
    elements_assembly : array
                        2d matrix which contains the ordering of the element based on the degree and knot span.
                        EX:[[0,1,2,3,4,5,6,7,8,9,10,11]]
    
    Test cases
    ----------
    test command -
                    pytest test_geometry.py::test__single_element_Assembly

                    pytest test_geometry.py::test__element_Assembly_C0_continuity

                    pytest test_geometry.py::test__single_element_Assembly_C1_continuity_along_x

    '''
    #Inititation of dimentions for 2D,3D matrices
    elements_assembly=np.zeros(((nU*nV*nW),(xdegree+1)*(ydegree+1)*(zdegree+1)))
    elements_order=elementorder(nx,ny,nz)
    a=0 
    #looped over the length of the knot in xi,eta,neta direction : No-of elements
    for i in range(nW):
        for j in range(nV):
            for k in range(nU):
                c=0
                #looped over knot span in xi,eta,neta direction : No-of knots within each element 
                for l in range(len(knotconnectivityW[i,:])):
                    for m in range(len(knotconnectivityV[j,:])):
                        for n in range(len(knotconnectivityU[k,:])):
                            elements_assembly[a,c]=elements_order[knotconnectivityU[k,n],knotconnectivityW[i,l], knotconnectivityV[j,m]]  #element index are generated based on knot index and element order
                            c+=1
                a+=1         
    elements_assembly=elements_assembly.astype(int)
    return elements_assembly

