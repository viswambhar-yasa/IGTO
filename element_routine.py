#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#--------------------------------------------------------------------------#
#ELEMENT ROUTINE - Python file used to generate element routine used in FEM 
#--------------------------------------------------------------------------#


import numpy as np
from geometry import controlpointassembly,knot_connectivity,trilinear_der,bspline_basis,derbspline_basis,knot_index
import pytest

def gauss_quadrature(p,q,r):
    '''
    Generates gauss newton points and weight based on sample points.

    Parameters
    ----------
    p,q,r : int
            number of sample points along x,y and z direction

    Returns
    -------
    gausspoints : Array
                  A 2D array conatining the gauss points 
                  Ex- [[0 0 0]]
    weights : Array
              An array which conatins the weights of a set of weight function.
              Ex- [8]
    
    Test cases 
    ----------

    test command -
                    pytest test_element_routine.py::test__gauss_quadrature_1point_true

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
        r=r+1
    
    [xgauss_points,xweights]=np.polynomial.legendre.leggauss(p)
    [ygauss_points,yweights]=np.polynomial.legendre.leggauss(q)
    [zgauss_points,zweights]=np.polynomial.legendre.leggauss(r)
    gausspoints=np.zeros((p*q*r,3))
    weights=np.zeros(p*q*r)
    n=0
    #looping over the sample points 
    for i in range(r):
        for j in range(q):
            for k in range(p):
                gausspoints[n,:]=[xgauss_points[k],ygauss_points[j],zgauss_points[i]] #gauss points 
                weights[n]=xweights[k]*yweights[j]*zweights[i]                        #gauss weights
                n+=1
    return gausspoints,weights



def unittoparametric(gauss_point,span):
    '''
    Used to generate parameters to build jacobian which maps parametric space to unit space 

    Parameters
    ----------
    gauss_point : float
                    Gauss newton points obtained from gauss quadrature
    span : array
                    knot span - knot values of the element 

    Returns
    -------
    Scalar value 

    Test cases 
    ----------
    test command - 
                    pytest test_element_routine.py::test__unit_to_parametric_space_true
    '''
    return (np.dot(0.5,((span[1]-span[0])*gauss_point+(span[1]+span[0]))))   #based on equ. 5.3,5.4,5.5



def jacobian(Xi,Eta,Nta,Uspan,Vspan,Wspan,X,Y,Z,weights,xdegree,xknot_vector,ydegree,yknot_vector,zdegree,zknot_vector):
    '''
    This function creates jacobian which maps from global to parameter space.

    Parameters
    ----------
    Xi,Eta,Nta : float
                    Values obtained from unittoparametric function along xi,eta,neta direction. 

    Uspan,Vspan,Wspan : array
                         A 2d array containing the knot values of the element.
 
    X,Y,Z : array
                A 2d array containing the control point co-ordinates of x,y and z direction.

    weights : array
                An array which contains the weights of each control point (NURBS)

    xdegree,ydegree,zdegree : int
                                 degree of the knotvector along xi,eta,neta
                                     0-constant, 1-linear, 2-quadratic, 3-cubic
       
    xknot_vector,yknot_vector,zknot_vector : array
                                             knot vector along xi,eta,neta direction.

    Returns
    -------
    jacobian1 : array
                    A jacobian matrix which maps physical space to parameteric space

    det_jacobian1 : float
                    determinant of jacobian1 (physical space to parameteric space)

    det_jacobian2 : float   
                        determinant of jacobian 2  (parameteric space to unit space)    
    
    nurbs : array
                Trivariant function obtained from trilinear derivative function
    
    Test cases
    ---------
    Test commands (Sanity check)-
                                    pytest test_element_routine.py::test__Jacobian_patch_testing_rotation
                    
    '''
    #Derivatives of the equ. 5.3,5.4,5.5
    DXi_dxi=0.5*(Uspan[1]-Uspan[0])
    DEta_dEta=0.5*(Vspan[1]-Vspan[0])
    DNta_dNta=0.5*(Wspan[1]-Wspan[0])
    #determinant of Jacobian2 (parameteric space to unit space) is calculated.
    det_jacobian2=DXi_dxi*DEta_dEta*DNta_dNta   # based on equ. 5.6
    pts=np.array([X,Y,Z])  #physical co-odinates system
    #Trivariant function and its derivatives are obtained from Trivariant 
    [dR_dxi,dR_deta,dR_dnta,nurbs]=trilinear_der(Xi,Eta,Nta,weights,xdegree,xknot_vector,ydegree,yknot_vector,zdegree,zknot_vector)
    DR_DXi=np.array([dR_dxi,dR_deta,dR_dnta])
    #Jacobian 1 (physical space to parameteric space)
    jacobian1=np.zeros((3,3)) 
    #Based on equ.5.2. The jacobian1 is calculated
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
    det_jacobian1=np.linalg.det(jacobian1) #determinant of jacobian1

    return jacobian1,det_jacobian1,det_jacobian2,DR_DXi,nurbs



def strain_displacement(Xi,Eta,Neta,jacobian1,det_jacobian1,DR_DXi,nurbs,nn):
    '''
    B matrix a constitutive matrix which relates stresses and strains.

    Parameters
    ----------
    Xi,Eta,Neta : float
                    Knot values obtained from unittoparametric function along xi,eta,neta direction. 

    jacobian1 : array
                    Maps physical space to parametric space obtained from jacobian function.

    det_jacobian1 : float
                    Determinant of jacobian1.
    DR_DXi : array
                Dervatives of trivariant function w.r.t xi,eta,neta.
    
    NURBS : array
                Trivariant functions.


    Returns
    -------
    B : array
        Strain displacement matrix .

    R : array
        Basis function .
    '''

    inverse_jacobian1=np.linalg.inv(jacobian1)
    DR_DX=np.matmul(inverse_jacobian1,DR_DXi)  #based on equ. 5.8 
    #initialization of the dimensions for B and R matrix
    B=np.zeros((6,int(3*np.size(DR_DX,1))))
    R=np.zeros((3,int(3*np.size(DR_DX,1))))

    # B and R matrix are calculated equ. 4.12 and 4.8
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
        R[1,y]=nurbs[i]
        R[2,z]=nurbs[i]
    return B,R


def Compliance_matrix(E=100000,v=0.3):
    '''
    Generates a isotropic elastic compliance tensor

    Parameters
    ----------
    E : int
        Young's modulus.
    v : float
        Poission's ratio.

    Returns
    -------
    C : array
        compliance tensor.

    '''
    #initialization of the dimensions for compliance tensor
    C=np.zeros((6,6))
    #6x6 compliance tensor is generated as it is isotropic
    C[0,0]=C[1,1]=C[2,2]=(E/((1+v)*(1-2*v)))*(1-v)
    C[3,3]=C[4,4]=C[5,5]=E/((1+v)*1)
    C[0,1]=C[0,2]=C[1,0]=C[2,0]=C[1,2]=C[2,1]=(E/((1+v)*(1-2*v)))*v
    return C

def element_routine(X,Y,Z,weights,E,v,Uspan,Vspan,Wspan,xdegree,xknot_vector,ydegree,yknot_vector,zdegree,zknot_vector):
    '''
    Element routine builds element stiffness matrix using jacobian and B matrix over gauss quadrature scheme.

    Parameters
    ----------
    Xi,Eta,Nta : float
                    Values obtained from unittoparametric function along xi,eta,neta direction. 

    weights : array
                An array which contains the weights of each control point (NURBS)


    Uspan,Vspan,Wspan : array
                         A 2d array containing the knot values of the element.
 
    X,Y,Z : array
                A 2d array containing the control point co-ordinates of x,y and z direction.


    xdegree,ydegree,zdegree : int
                                 degree of the knotvector along xi,eta,neta
                                     0-constant, 1-linear, 2-quadratic, 3-cubic
       
    xknot_vector,yknot_vector,zknot_vector : array
                                             knot vector along xi,eta,neta direction.

    Returns
    -------
    Ke : array
                elementary stiffness matrix

    nurbs : array
                Trivariant basis functions

    R : array   
              Trivariant basis functions mapped according to equ. 4.8  obtained from jacobian
    
    B : array
            strain displacement matrix  obtained from jacobian
    
    Test cases
    ---------
    Test commands (Sanity check)-
                                    pytest test_element_routine.py::test__stiffness_matrix_singularity
    '''
    #Gauss quadrature is obtained based on knot degree along xi,eta,neta
    Gauss_points,gauss_weights=gauss_quadrature(xdegree,ydegree,zdegree)
    C=Compliance_matrix(E,v)  #compliance tensor is obtained
    #initialization of dimensions for element stiffness matrix 
    Ke=np.zeros((3*((xdegree+1)*(ydegree+1)*(zdegree+1)),3*((xdegree+1)*(ydegree+1)*(zdegree+1))))
    nn=((xdegree+1)*(ydegree+1)*(zdegree+1)) #number of control points in each element
    #looped over gauss quadature
    for i in range(len(gauss_weights)):
        gp=Gauss_points[i,:]
        Wt=gauss_weights[i]
        Xi=unittoparametric(gp[0],Uspan) #based on equ. 5.3
        Eta=unittoparametric(gp[1],Vspan) #based on equ 5.4 
        Neta=unittoparametric(gp[2],Wspan) #based on equ. 5.5
        #jacobian are obtained 
        [jacobian1,det_jacobian1,det_jacobian2,DR_DXi,nurbs]=jacobian(Xi,Eta,Neta,Uspan,Vspan,Wspan,X,Y,Z,weights,xdegree,xknot_vector,ydegree,yknot_vector,zdegree,zknot_vector)
        #strain-displacement matrix is obtained
        B,R=strain_displacement(Xi,Eta,Neta,jacobian1,det_jacobian1,DR_DXi,nurbs,nn)
        temp=np.transpose(B)@C
        #element stiffness matrix is built using strain-displacement and jacobians
        Ke+=temp@B*(det_jacobian1*det_jacobian2*Wt)  #based on equ. 5.7
    return Ke,nurbs,R,B

def stress_strain_element_routine(X,Y,Z,weights,displacements,E,v,Uspan,Vspan,Wspan,xdegree,xknot_vector,ydegree,yknot_vector,zdegree,zknot_vector):
    '''
    This function builds stress and strains at gauss quadrature based on the displacement 

    Parameters
    ----------
    X,Y,Z : array
                A 2d array containing the control point co-ordinates of x,y and z direction. 

    weights : array
                An array which contains the weights of each control point (NURBS)


    Uspan,Vspan,Wspan : array
                         A 2d array containing the knot values of the element.
 
    E :int
            Young's modulus   

    v : float
            poission's ratio

    xdegree,ydegree,zdegree : int
                                 degree of the knotvector along xi,eta,neta
                                     0-constant, 1-linear, 2-quadratic, 3-cubic
       
    xknot_vector,yknot_vector,zknot_vector : array
                                             knot vector along xi,eta,neta direction.

    Returns
    -------
    gp_strain : array
                element strains at gauss quadrature

    gp_stress : array
                element stress at gauss quadrature

    R : array   
              Trivariant basis functions mapped according to equ. 4.8  obtained from jacobian
    
    B : array
            strain displacement matrix obtained from jacobian

    '''
    nn=((xdegree+1)*(ydegree+1)*(zdegree+1)) #number of control point in each element
    C=Compliance_matrix(E,v)  #Compliance tensor
    #initialization of dimensions for gauss point stresses and strain matrix
    gp_strain=np.zeros((nn,6))
    gp_stress=np.zeros((nn,6))
    gp_disp=np.zeros((nn,3))
    gp=0
    pts=np.array([X,Y,Z]) #co-ordinates of control points

    u_index=knot_index(xdegree,Uspan[0],xknot_vector)
    v_index=knot_index(ydegree,Vspan[0],yknot_vector)
    w_index=knot_index(zdegree,Wspan[0],zknot_vector)
    #looped over the span of knots (no of control point each element along xi,eta,neta)
    for i in range(len(Wspan)):
        for j in range(len(Vspan)):
            for k in range(len(Uspan)):
                Xi=Uspan[k]
                Eta=Vspan[j]
                Neta=Wspan[i]
                # Basis function along xi,eta,neta
                Nx=bspline_basis(u_index,xdegree,Xi,xknot_vector)
                Ny=bspline_basis(v_index,ydegree,Eta,yknot_vector)
                Nz=bspline_basis(w_index,zdegree,Neta,zknot_vector)

                #derivatives of basis function along xi,eta,neta
                DNx=derbspline_basis(u_index,xdegree,Xi,xknot_vector)
                DNy=derbspline_basis(v_index,ydegree,Eta,yknot_vector)
                DNz=derbspline_basis(w_index,zdegree,Neta,zknot_vector)

                W=0
                W_dx=0
                W_dy=0
                W_dz=0
                p=0
                windex=0
                #calculating weights to use in trivariant function 
                for k in range(zdegree+1):
                    for j in range(ydegree+1):
                        for i in range(xdegree+1):
                            W+=Nx[i]*Ny[j]*Nz[k]*weights[windex]     #Based on equ. 3.9
                            W_dx+=DNx[i]*Ny[j]*Nz[k]*weights[windex] #Based on equ. 3.10
                            W_dy+=Nx[i]*DNy[j]*Nz[k]*weights[windex] #Based on equ. 3.11
                            W_dz+=Nx[i]*Ny[j]*DNz[k]*weights[windex] #Based on equ. 3.12
                            windex+=1

                dR_dx=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
                dR_dy=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
                dR_dz=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))
                R=np.zeros((xdegree+1)*(ydegree+1)*(zdegree+1))

                w=0
                windex=0
                #building trivariant basis function w.r.t xi,eta,neta
                for k in range(zdegree+1):
                    for j in range(ydegree+1):
                        for i in range(xdegree+1):
                            w=weights[windex]/(W*W)
                            R[p]=(Nx[i]*Ny[j]*Nz[k]*W*w)                                #Based on equ. 3.5
                            dR_dx[p]=(DNx[i]*Ny[j]*Nz[k]*W-Nx[i]*Ny[j]*Nz[k]*W_dx)*w    #Based on equ. 3.6
                            dR_dy[p]=(Nx[i]*DNy[j]*Nz[k]*W-Nx[i]*Ny[j]*Nz[k]*W_dy)*w    #Based on equ. 3.7
                            dR_dz[p]=(Nx[i]*Ny[j]*DNz[k]*W-Nx[i]*Ny[j]*Nz[k]*W_dz)*w    #Based on equ. 3.8
                            p=p+1
                            windex+=1
                #derivative of trivariant basis function w.r.t x,y,z
                DR_DXi=np.array([dR_dx,dR_dy,dR_dz])
                jacobian=pts@ np.transpose(DR_DXi)
                det_jacobian=np.linalg.det(jacobian)
                B,R=strain_displacement(Xi,Eta,Neta,jacobian,det_jacobian,DR_DXi,R,nn)
                # Stresses and strains at gauss points are calculated.
                gp_strain[gp,:]=np.matmul(B,displacements)
                gp_stress[gp,:]=np.matmul(C,gp_strain[gp,:])
                gp_disp[gp,:]=np.matmul(R,displacements.T)
                gp+=1
    M=stress_interploation_matrix()   
    #nodal_stress=gp_stress
    #nodal_strain=gp_strain 
    return gp_strain,gp_stress,gp_disp

def stress_interploation_matrix():
    '''
    Maps stress and strains calculated at gauss points to nodal points
    '''
    a=(5+3*(3**(0.5)))/4
    b=(1+3**(0.5))/4
    c=(3**(0.5)-1)/4
    d=(5-3*(3**(0.5)))/4
    M=np.array([[a,b,b,c,b,c,c,d],
                [b,c,c,d,a,b,b,c],
                [c,d,b,c,b,c,a,b],
                [b,c,a,b,c,d,b,c],
                [b,a,c,b,c,b,d,c],
                [c,b,d,c,b,a,c,b],
                [d,c,c,b,c,b,b,a],
                [c,b,b,a,d,c,c,b]])
    return M


def assemble(K_G,K_E,elindices,ncp,K_disp=False):
    '''
    Builds Global stiffness matrix from element stiffness matrix based on element indicies(similar to Assignment matrix).

    Parameters
    ----------
    K_G : array
        Global stiffness matrix.
    K_E : array
        Element stiffness matrix obtained from element routine.
    elindices : array
            An array containg the element indices to map element stiffness into global stiffness.
    ncp : int
        Number of control point in each element.

    Returns
    -------
    K_G : array
        Global stifness matrix.

    '''
    
    elindices=np.array(elindices)
    dof=3
    gindices=np.sort(np.concatenate((elindices*dof,dof*elindices+1,dof*elindices+2)))
    #looped over element indices
    for i,ii in enumerate(gindices):
        for j,jj in enumerate(gindices):
            K_G[ii,jj]=K_G[ii,jj]+K_E[i,j]  #based on equ. 4.17
    if K_disp:        
        print('Building Global Stiffness matrix \n')
    return K_G


def apply_BC(F_E,fixed_dof,load_dof,P,option=0,abc_disp=False):
    '''
    This function is used to apply Boundary conditions.
    IGA requires additional method to apply boundary condition like traction, so seperate function is generated for future scope. 

    Parameters
    ----------
    F_E : array
        An empty external force vector.

    fixed_dof : array
        An array which contains element indicies which are fixed 

    load_dof : array
        An array which contains element indicies where load is appled.
    
    P : float
        load or force.

    Returns
    -------
    reduced_F_E : array
                    External force vector after Boundary conditions are applied.

    '''
    if option==4:
        F_E[load_dof[0]]=P
        F_E[load_dof[1]]=-P
    else:
        F_E[load_dof]=P
    reduced_F_E=np.delete(F_E,fixed_dof,axis=0)
    if abc_disp:
        print('Applying Boundary Conditions \n')
    return reduced_F_E
