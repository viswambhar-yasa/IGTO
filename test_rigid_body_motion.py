
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
    Using eigen value test , we find the rigid body motion 
    The number of eigen value which are equal to zero represent rigid body motions
    Ex- 3 (for 2d element) 
    2 translation and 1 rotation

    '''
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

    ncp = N * P * Q
    dof = 3
    dofcp = ncp * dof
    nel = nU * nV * nW

    K_G = np.zeros((dofcp, dofcp))
    F_E = np.zeros(dofcp)
    K_disp = True
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
    ncp = N * P * Q
    dof = 3
    dofcp = ncp * dof
    nel = nU * nV * nW
    K_G = np.zeros((dofcp, dofcp))
    F_E = np.zeros(dofcp)
    U = np.zeros(dofcp)
    K_disp = True
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

    F=np.matmul(K_G,U)
    disp_dof=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])
    disp_index=np.sort((disp_dof*3))
    #assigning outer node with displacement
    U[disp_index]=0.5
    fixed_indicies=np.sort(np.concatenate((3*disp_dof+1,disp_dof*3+2)))
    #assigning inner node with no displacement
    U[13*3]=0
    #calculating the residual force due to displacement along x axis
    residual_F=np.matmul(K_G,U)
    '''
    load_dof=np.array([14])
    fixed_indicies=np.sort(load_dof*3+1)
    F[fixed_indicies]=0
    '''
    #calculating the displacement along in the inner node 
    U[13*3]=(F[13*3]-residual_F[13*3])/(K_G[13*3,13*3])
    assert (all(U[disp_index]==0.5)) is True



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
    CONTROL_POINTS[13,:]=np.array([ 0.6,0.4,0.3,1.,13.])
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

