#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#--------------------------------------------------------------------------------#
#A python test file to run rigid body motion (validating global stiffness matrix)
# command to run all test cases
# pytest test_rigid_body_motion.py
# -------------------------------------------------------------------------------# 
from Preprocessing import *
from geometry import knot_connectivity, controlpointassembly
from element_routine import assemble, element_routine, apply_BC,Compliance_matrix ,stress_strain_element_routine,gauss_quadrature
from boundary_conditions import BC_Switcher
from analytical_solution import exact_displacements
import math as m
import time 
import pytest
import random


def test__rigid_body_constrains():
    '''
    PATCH TESTING
    Aim : Also known as eigen value test ,We check rigid body motion for constrained and unconstrained global stiffness matrix
    
    Expected Output : The number of zero eigenvalues from global stiffness matrix represent the rigid body motion of the structure.
                    unconstrained - 6 zeros eigenvalues (3-rotation, 3-translation)
                    constrained - 0 zero eigenvalues
                    Ex- 3 (for 2d element) 
                    2 translation and 1 rotation
    
    Test command : pytest test_rigid_body_motion.py::test_rigid_body_constrains

    Remarks : Test case passed.

    '''
    #input parameter are initialized
    length=1
    height=1
    width=1
    nx=5
    ny=5
    nz=5

    XI_DEGREE=1
    ETA_DEGREE=1
    NETA_DEGREE=1

    Youngs_modulus=100000
    poission_ratio=0.3

    N = nx
    P = ny
    Q = nz
    option=3
    load=-200
    #knot vector and knot span are obtained from input parameter     
    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS = C.crtpts_coordinates()
    # print(CONTROL_POINTS)
    WEIGHTS = CONTROL_POINTS[:, -1]

    XI_KNOTVECTOR = C.xi_knotvector()
    ETA_KNOTVECTOR = C.eta_knotvector()
    NETA_KNOTVECTOR = C.neta_knotvector()
    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()

    ncp = N * P * Q # no of control points
    dof = 3
    dofcp = ncp * dof # nof of nodal degree of freedom
    nel = nU * nV * nW  
    
    K_G = np.zeros((dofcp, dofcp))
    F_E = np.zeros(dofcp)
    K_disp = True
    #Control point assembly is build based on knot span and knot connectivity
    element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                            ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    span_index = knot_connectivity(N, P, Q, XI_KNOTCONNECTIVITY, ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    #looped over number of elements
    for i in range(0, nel):
        el_in = element_indicies[i, :]
        sp_in = span_index[i, :]
        #co-ordinate of control points in physical space
        X = CONTROL_POINTS[el_in, 0]
        Y = CONTROL_POINTS[el_in, 1]
        Z = CONTROL_POINTS[el_in, 2]
        weights = CONTROL_POINTS[el_in, 3]
        Uspan = XI_SPAN[sp_in[0], :]
        Vspan = ETA_SPAN[sp_in[1], :]
        Wspan = NETA_SPAN[sp_in[2], :]
        #Element stiffness matrix is obtained
        K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                        XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)
        #Global stiffness matrix is assembled from element stiffness matrix using element indices from control point assembly for mapping.
        K_G = assemble(K_G, K_E, el_in, ncp, K_disp)
        K_disp = False
    #eigen values are calculated for global stiffness matrix
    eigen_values,eigen_vectors= np.linalg.eig(K_G)
    eigen_values=np.round(eigen_values,8)
    unconstrained_rigid_body_modes=len(eigen_values[eigen_values==0])
    #For 3d unconstrained structure
    #3 translation and 3 rotation should be present
    # total of 6 zero eigen values should be present for stiffness matrix
    bc_disp = False
    abc_disp = True
    BC = BC_Switcher(CONTROL_POINTS, length, height, width, bc_disp)
    fixed_dof, load_dof, fixed_pts, load_pts = BC.indirect(option)
    reduced_k=np.delete(np.delete(K_G, fixed_dof, 0),fixed_dof , 1)
    reduced_F = apply_BC(F_E, fixed_dof, load_dof, load, abc_disp)
    eigen_values,eigen_vectors= np.linalg.eig(reduced_k)
    eigen_values=np.round(eigen_values,8)
    #print(eigen_values)
    constrained_rigid_body_modes=len(eigen_values[eigen_values==0])
    #Too few zero eigenvalues is an indication of an element lacking the capability of rigid body motion without strain.
    #Translation along y axis Cantilever beam
    assert (unconstrained_rigid_body_modes==6) and (constrained_rigid_body_modes==0) is True


def test__rigid_body_translation_along_x():
    '''
    PATCH TESTING
    We take a 3 x 3 x 3 structure with interior nodes, we apply displacement boundary condition on all exterior nodes and calculate the residual force. From the residual forces , the interior node displacements are calculated .
    
    Aim : Rigid body translation along x direction-global stiffness matrix should be able to represent rigid body translation.

    Expected Output : The displacement of the external and internal nodes should be same.
    
    Test command : pytest test_rigid_body_motion.py::test_rigid_body_translation_along _x
    
    Remarks : Test case passed.
    '''
    #input parameter are initialized
    length=1
    height=1
    width=1
    nx=3
    ny=3
    nz=3
    XI_DEGREE=1
    ETA_DEGREE=1
    NETA_DEGREE=1
    Youngs_modulus=100000
    poission_ratio=0.3
    N = nx
    P = ny
    Q = nz
    #knot vector and knot span are obtained from input parameter 
    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)
    CONTROL_POINTS = C.crtpts_coordinates()
    CONTROL_POINTS[13,:]=np.array([ 0.6,0.4,0.3,1.,13.])
    WEIGHTS = CONTROL_POINTS[:, -1]
    XI_KNOTVECTOR = C.xi_knotvector()
    ETA_KNOTVECTOR = C.eta_knotvector()
    NETA_KNOTVECTOR = C.neta_knotvector()
    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()
    ncp = N * P * Q # no of control points
    dof = 3
    dofcp = ncp * dof # nof of nodal degree of freedom
    nel = nU * nV * nW #no of elements
    
    
    K_G = np.zeros((dofcp, dofcp))
    F_E = np.zeros(dofcp)
    U = np.zeros(dofcp)
    K_disp = True
    #Control point assembly is build based on knot span and knot connectivity
    element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                            ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    span_index = knot_connectivity(N, P, Q, XI_KNOTCONNECTIVITY, ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    #looped over number of elements
    for i in range(0, nel):
        el_in = element_indicies[i, :]
        sp_in = span_index[i, :]
        #co-ordinate of control points in physical space
        X = CONTROL_POINTS[el_in, 0]
        Y = CONTROL_POINTS[el_in, 1]
        Z = CONTROL_POINTS[el_in, 2]
        weights = CONTROL_POINTS[el_in, 3]
        Uspan = XI_SPAN[sp_in[0], :]
        Vspan = ETA_SPAN[sp_in[1], :]
        Wspan = NETA_SPAN[sp_in[2], :]
        #Element stiffness matrix is obtained
        K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                        XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)
        ##Global stiffness matrix is assembled from element stiffness matrix using element indices from control point assembly for mapping
        K_G = assemble(K_G, K_E, el_in, ncp, K_disp)
        K_disp = False
    #Resdiual force due to displacement of external nodes is calculated
    F=np.matmul(K_G,U)
    #displacement based boundary condition are formulated
    disp_dof=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])
    disp_index=np.sort((disp_dof*3))
    #assigning outer node with displacement
    U[disp_index]=0.5
    fixed_indicies=np.sort(np.concatenate((3*disp_dof+1,disp_dof*3+2)))
    #assigning inner node with no displacement
    U[13*3]=0
    #calculating the residual force due to displacement along x axis
    residual_F=np.matmul(K_G,U)
    
    #calculating the displacement along in the inner node 
    U[13*3]=(F[13*3]-residual_F[13*3])/(K_G[13*3,13*3])
    #Inner and outer node displacement should match
    assert (all(U[disp_index]==0.5)) is True



def test__rigid_body_rotation():
    '''
    we use Euler Bunge transformation matrix
    The external nodes are rotation and internal nodes displacement are calculated.
    Aim : Rigid body rotation - global stiffness matrix should be able to represent rigid body rotation .
    
    Expected Output : The displacement of the internal node should match the coordinates obtained after performing rotation using equ.(6.7) to reference system.
    
    Test command : pytest test_rigid_body_motion.py::test_rigid_body_rotation
    
    Remarks : Test case passed.
    '''
    #input parameters are initialized
    length=1
    height=1
    width=1
    nx=3
    ny=3
    nz=3
    Youngs_modulus=100000
    poission_ratio=0.3
    XI_DEGREE = 1
    ETA_DEGREE = 1
    NETA_DEGREE = 1
    N = nx
    P = ny
    Q = nz
    #knot vector and knot span are obtained from input parameter
    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)
    
    CONTROL_POINTS=C.crtpts_coordinates()
    #generated an irregural control points
    CONTROL_POINTS[13,:]=np.array([ 0.6,0.4,0.3,1.,13.])
    XI_KNOTVECTOR = C.xi_knotvector()
    ETA_KNOTVECTOR = C.eta_knotvector()
    NETA_KNOTVECTOR = C.neta_knotvector()
    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()
    ncp = N * P * Q # no of control points
    dof = 3
    dofcp = ncp * dof # no of nodal degree of freedom
    nel = nU * nV * nW 
    #initializing dimension of global stiffness matrix and external force vector
    K_G = np.zeros((dofcp, dofcp))
    F_E = np.zeros(dofcp)
    U = np.zeros(dofcp)
    K_disp = False
    #Control point assembly is build based on knot span and knot connectivity
    element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                            ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    span_index = knot_connectivity(N, P, Q, XI_KNOTCONNECTIVITY, ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    #looped over number of elements
    for i in range(0, nel):
        el_in = element_indicies[i, :]
        sp_in = span_index[i, :]
        #co-ordinate of control points in physical space
        X = CONTROL_POINTS[el_in, 0]
        Y = CONTROL_POINTS[el_in, 1]
        Z = CONTROL_POINTS[el_in, 2]
        weights = CONTROL_POINTS[el_in, 3]
        Uspan = XI_SPAN[sp_in[0], :]
        Vspan = ETA_SPAN[sp_in[1], :]
        Wspan = NETA_SPAN[sp_in[2], :]
        #Element stiffness matrix is obtained
        K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                        XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)
         ##Global stiffness matrix is assembled from element stiffness matrix using element indices from control point assembly for mapping
        K_G = assemble(K_G, K_E, el_in, ncp, K_disp)
        K_disp = False
    ctrl_pts=CONTROL_POINTS[:,:-2]

    #print(ctrl_pts)           
    a1=np.radians(60)
    a2=np.radians(0)
    a3=np.radians(30)
    #rotation matrix
    Q1=np.array([[m.cos(a1),m.sin(a1) ,0],[-m.sin(a1),m.cos(a1) ,0],[0,0,1]])
    Q2=np.array([[m.cos(a2),0,-m.sin(a2)],[0,1,0],[m.sin(a2),0,m.cos(a2)]])
    Q3=np.array([[m.cos(a3),m.sin(a3) ,0],[-m.sin(a3),m.cos(a3) ,0],[0,0,1]])
    #Transformation matrix using Euler Bunge formulation
    rotation_matrix=Q1@Q3@Q1 #(Z-X-Z)
    New_control_points=np.matmul(ctrl_pts,rotation_matrix)
    for i in range(0,len(CONTROL_POINTS[:,0])):
        if i!=13:#assigning displacements of outer nodes
            disp_index=np.array([i*3,i*3+1,i*3+2])
            U[disp_index]=np.matmul(ctrl_pts[i,:],rotation_matrix)
        else:
            disp_index=np.array([i*3,i*3+1,i*3+2])
            U[disp_index]=np.matmul(np.array([0,0,0]),rotation_matrix)
    #residual forces are calculated after the exterior nodes are rotated        
    residual_F=np.matmul(K_G,U)
    #displacement based on Boundary conditon are formulated.
    fixed_dof=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26])
    fixed_indicies=np.sort(np.concatenate((3*fixed_dof,3*fixed_dof+1,fixed_dof*3+2)))
    F=F_E-residual_F
    reduced_k=np.delete(np.delete(K_G, fixed_indicies, 0),fixed_indicies , 1)
    reduced_F = np.delete(F,fixed_indicies)
    #calculating the displacements of inner node
    U_disp=np.linalg.solve(reduced_k,reduced_F)
    #print(U_disp)
    U_inner_node_exact=New_control_points[13]
    #inner node displacement calculated by FEM should match with the inner node obtained after rotating using Euler Bunge transformation matrix
    output=np.all(abs(U_disp-U_inner_node_exact)<1e-5)
    o=0
    if output == True:
        o=1
    assert (o==1) is True


def test__rigid_body_rotation():
    length=1
    height=1
    width=1
    nx=3
    ny=3
    nz=3
    Youngs_modulus=100000
    poission_ratio=0.3
    XI_DEGREE = 1
    ETA_DEGREE = 1
    NETA_DEGREE = 1
    N = nx
    P = ny
    Q = nz
    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS=C.crtpts_coordinates()
    XI_KNOTVECTOR = C.xi_knotvector()
    ETA_KNOTVECTOR = C.eta_knotvector()
    NETA_KNOTVECTOR = C.neta_knotvector()
    XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
    ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
    NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()
    ncp = N * P * Q
    dof = 3
    dofcp = ncp * dof
    nel = nU * nV * nW
    K_G = np.zeros((dofcp, dofcp))
    F_E = np.zeros(dofcp)
    U = np.zeros(dofcp)
    K_disp = False
    element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                            ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    span_index = knot_connectivity(N, P, Q, XI_KNOTCONNECTIVITY, ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
    for i in range(0, nel):
        el_in = element_indicies[i, :]
        sp_in = span_index[i, :]
        X = CONTROL_POINTS[el_in, 0]
        Y = CONTROL_POINTS[el_in, 1]
        Z = CONTROL_POINTS[el_in, 2]
        weights = CONTROL_POINTS[el_in, 3]
        Uspan = XI_SPAN[sp_in[0], :]
        Vspan = ETA_SPAN[sp_in[1], :]
        Wspan = NETA_SPAN[sp_in[2], :]
        K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                        XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)
        K_G = assemble(K_G, K_E, el_in, ncp, K_disp)
        K_disp = False
    ctrl_pts=CONTROL_POINTS[:,:-2]

    #print(ctrl_pts)           
    a1=np.radians(60)
    a2=np.radians(0)
    a3=np.radians(30)
    Q1=np.array([[m.cos(a1),m.sin(a1) ,0],[-m.sin(a1),m.cos(a1) ,0],[0,0,1]])
    Q2=np.array([[m.cos(a2),0,-m.sin(a2)],[0,1,0],[m.sin(a2),0,m.cos(a2)]])
    Q3=np.array([[m.cos(a3),m.sin(a3) ,0],[-m.sin(a3),m.cos(a3) ,0],[0,0,1]])
    rotation_matrix=Q1@Q3@Q1 #(Z-X-Z)
    New_control_points=np.matmul(ctrl_pts,rotation_matrix)
    for i in range(0,len(CONTROL_POINTS[:,0])):
        if i!=13:#assigning displacements of outer nodes
            disp_index=np.array([i*3,i*3+1,i*3+2])
            U[disp_index]=np.matmul(ctrl_pts[i,:],rotation_matrix)
        else:
            disp_index=np.array([i*3,i*3+1,i*3+2])
            U[disp_index]=np.matmul(np.array([0,0,0]),rotation_matrix)
    residual_F=np.matmul(K_G,U)
    #calculating the displacements of inner node
    fixed_dof=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26])
    fixed_indicies=np.sort(np.concatenate((3*fixed_dof,3*fixed_dof+1,fixed_dof*3+2)))
    F=F_E-residual_F
    reduced_k=np.delete(np.delete(K_G, fixed_indicies, 0),fixed_indicies , 1)
    reduced_F = np.delete(F,fixed_indicies)
    
    U_disp=np.linalg.solve(reduced_k,reduced_F)
    #print(U_disp)
    U_inner_node_exact=New_control_points[13]
    output=np.all(abs(U_disp-U_inner_node_exact)<1e-5)
    o=0
    if output == True:
        o=1
    assert (o==1) is True


tol=1e-2
length=1
height=1
width=1
nx=2
ny=2
nz=2
Youngs_modulus=100000
poission_ratio=0.3
XI_DEGREE = 1
ETA_DEGREE = 1
NETA_DEGREE = 1
N = nx
P = ny
Q = nz
C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)
CONTROL_POINTS=C.crtpts_coordinates()
XI_KNOTVECTOR = C.xi_knotvector()
ETA_KNOTVECTOR = C.eta_knotvector()
NETA_KNOTVECTOR = C.neta_knotvector()
XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()
ncp = N * P * Q
dof = 3
dofcp = ncp * dof
nel = nU * nV * nW
K_G = np.zeros((dofcp, dofcp))
K_G2 = np.zeros((dofcp, dofcp))
F_E = np.zeros(dofcp)
U = np.zeros(dofcp)
K_disp = False
element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                        ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
span_index = knot_connectivity(N, P, Q, XI_KNOTCONNECTIVITY, ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)

for i in range(0, nel):
    el_in = element_indicies[i, :]
    sp_in = span_index[i, :]
    X = CONTROL_POINTS[el_in, 0]
    Y = CONTROL_POINTS[el_in, 1]
    Z = CONTROL_POINTS[el_in, 2]
    weights = CONTROL_POINTS[el_in, 3]
    Uspan = XI_SPAN[sp_in[0], :]
    Vspan = ETA_SPAN[sp_in[1], :]
    Wspan = NETA_SPAN[sp_in[2], :]
    K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                    XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)
    K_G = assemble(K_G, K_E, el_in, ncp, K_disp)
    K_disp = False

print(CONTROL_POINTS)
fixed_nodes=[0,1,2,4,6,12,13]
load_index=np.array([1,3,5,7])
F_E[3*load_index]=1
#F_residual=np.matmul(K_G,U)

reduced_k=np.delete(np.delete(K_G, fixed_nodes, 0),fixed_nodes , 1)
reduced_F = np.delete((F_E),fixed_nodes)

U_disp=np.linalg.solve(reduced_k,reduced_F)
for j in fixed_nodes:
    U_disp = np.insert(U_disp, j, 0)
print(U_disp.reshape(8,3))  

fixed_nodes=[0,1,2,4,6,12,13]
load_index=np.array([1,3,5,7])
F_E[3*load_index]=1/4
angle=15           
a1=np.radians(15)
a2=np.radians(0)
a3=np.radians(0)
Q1=np.array([[m.cos(a1),m.sin(a1) ,0],[-m.sin(a1),m.cos(a1) ,0],[0,0,1]])
Q2=np.array([[m.cos(a2),0,-m.sin(a2)],[0,1,0],[m.sin(a2),0,m.cos(a2)]])
Q3=np.array([[m.cos(a3),m.sin(a3) ,0],[-m.sin(a3),m.cos(a3) ,0],[0,0,1]])
#Z-X-Z
rotation_matrix=np.matmul(Q1,np.matmul(Q2,Q3))
print(rotation_matrix,CONTROL_POINTS[:,:-2])
rotated_CONTROL_POINTS=np.matmul(CONTROL_POINTS[:,:-2],rotation_matrix)
for i in range(0, nel):
    el_in = element_indicies[i, :]
    sp_in = span_index[i, :]
    X = rotated_CONTROL_POINTS[el_in, 0]
    Y = rotated_CONTROL_POINTS[el_in, 1]
    Z = rotated_CONTROL_POINTS[el_in, 2]
    weights = CONTROL_POINTS[el_in, 3]
    Uspan = XI_SPAN[sp_in[0], :]
    Vspan = ETA_SPAN[sp_in[1], :]
    Wspan = NETA_SPAN[sp_in[2], :]
    K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                    XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)
    K_G2 = assemble(K_G, K_E, el_in, ncp, K_disp)
    K_disp = False
#rotation_U=np.matmul(U.reshape(8,3),rotation_matrix).flatten()
#print(rotation_U)
#rotated_F_residual=np.matmul(K_G,rotation_U)
#F_residual=np.matmul(K_G,U)
rotated_F_E=np.matmul(F_E.reshape(8,3),rotation_matrix).flatten()
print(rotated_F_E)
#print(F_E)
reduced_k=np.delete(np.delete(K_G2, fixed_nodes, 0),fixed_nodes , 1)
reduced_F = np.delete((rotated_F_E),fixed_nodes)

U_disp=np.linalg.solve(reduced_k,reduced_F)
for j in fixed_nodes:
    U_disp = np.insert(U_disp, j, 0)
print(U_disp.reshape(8,3)) 