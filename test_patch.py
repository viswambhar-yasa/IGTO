from test_preprocessing import *
from geometry import knot_connectivity, controlpointassembly
from element_routine import assemble, element_routine, apply_BC,Compliance_matrix ,stress_strain_element_routine,gauss_quadrature
from boundary_conditions import BC_Switcher
from analytical_solution import exact_displacements
import math as m
import time 
import pytest
import random


def test__patch_constant_stress():
    '''
    for tensile loading,the stress generated within the element should be same at each gauss point 2x2x2 (8 gauss points)
    Simple pacth test

    Constant stress and strain in the element 

    '''
    length=1
    height=1
    width=0.1
    nx=2
    ny=2
    nz=2
    Youngs_modulus=150000
    poission_ratio=0.35
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
    fixed_nodes=[0,1,2,4,6,12,13]
    load_index=np.array([1,3,5,7])
    U[3*load_index]=1
    print(U)
    F_residual=np.matmul(K_G,U)
    print(F_residual)

    reduced_k=np.delete(np.delete(K_G, fixed_nodes, 0),fixed_nodes , 1)
    reduced_F = np.delete((F_residual),fixed_nodes)

    U_disp=np.linalg.solve(reduced_k,reduced_F)

    for j in fixed_nodes:
         U_disp = np.insert(U_disp, j, 0)
         
    nn=(XI_DEGREE+1)*(ETA_DEGREE+1)*(NETA_DEGREE+1)
    gauss_points_strains=np.zeros((nel,nn,6))
    gauss_point_stresses=np.zeros((nel,nn,6))
    gauss_point_disp=np.zeros((nel,nn,3))
    C=Compliance_matrix(Youngs_modulus,poission_ratio)
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
        u_index=np.sort(np.concatenate((el_in*dof,dof*el_in+1,el_in*dof+2)))
        U_el=U_disp[u_index].T
        gauss_points_strains[i,:,:],gauss_point_stresses[i,:,:],gauss_point_disp[i,:,:] = stress_strain_element_routine(X, Y, Z, weights,U_el, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                        XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)                              
    stresses=gauss_point_stresses[0,:,:]
    sigma_x=stresses[:,0]
    sigma_y=stresses[:,1]
    sigma_z=stresses[:,2]
    constant_sigma=(np.all(abs(sigma_x-sigma_x[0])<1e-5) and np.all(abs(sigma_y-sigma_y[0])<1e-5) and np.all(abs(sigma_z-sigma_z[0])<1e-5))
    o=0
    if constant_sigma :
        o=1
    assert (o==1) is True


def test__patch_analytical_solution():
    tol=1e-2
    length=48
    height=12
    width=1
    option=3
    nx=27
    ny=13
    nz=3
    load=-20

    Youngs_modulus=100000
    poission_ratio=0.3

    XI_DEGREE = 1
    ETA_DEGREE = 1
    NETA_DEGREE = 1

    N = nx
    P = ny
    Q = nz

    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS = C.crtpts_coordinates()
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

    bc_disp = False
    abc_disp = False
    BC = BC_Switcher(CONTROL_POINTS, length, height, width, bc_disp)
    fixed_dof, load_dof, fixed_pts, load_pts = BC.indirect(option)

    reduced_k=np.delete(np.delete(K_G, fixed_dof, 0),fixed_dof , 1)
    reduced_F = apply_BC(F_E, fixed_dof, load_dof, load, abc_disp)


    U = np.linalg.solve(reduced_k, reduced_F)

    for j in fixed_dof:
        U = np.insert(U, j, 0)

    F_E[load_dof] = load

    U_disp=np.concatenate((CONTROL_POINTS[:,:-2], U.reshape(len(CONTROL_POINTS), 3)),axis=1)
    analy_DISP_Y=[]
    Disp_y=[]
    analytical_solution = np.array(exact_displacements(load, length, width, height, Youngs_modulus, poission_ratio, nx, ny,nz))
    #deflection at  point load
    for i in range(len(analytical_solution[:,0])):
        if analytical_solution[i,2]==width/2 and analytical_solution[i,0]==length and analytical_solution[i,1]==height:
            analy_DISP_Y.append(analytical_solution[i,4])
        if U_disp[i,2]==width/2 and U_disp[i,0]==length and U_disp[i,1]==height:   
            Disp_y.append(U_disp[i,4])
    #deflection at  center    
    for i in range(len(analytical_solution[:,0])):
        if analytical_solution[i,2]==0 and analytical_solution[i,0]==length/2 and analytical_solution[i,1]==height/2:
            analy_DISP_Y.append(analytical_solution[i,4])
        if U_disp[i,2]==0 and U_disp[i,0]==length/2 and U_disp[i,1]==height/2:   
            Disp_y.append(U_disp[i,4])
    #deflection at bottom end
    for i in range(len(analytical_solution[:,0])):
        if analytical_solution[i,2]==0 and analytical_solution[i,0]==length and analytical_solution[i,1]==0:
            analy_DISP_Y.append(analytical_solution[i,4])
        if U_disp[i,2]==0 and U_disp[i,0]==length and U_disp[i,1]==0:   
            Disp_y.append(U_disp[i,4])
    expected_output=np.array(analy_DISP_Y)
    output=np.array(Disp_y)
    er=np.all(abs(expected_output-output)<tol)
    o=0
    if er:
        o=1
    assert (o==1) is True

def test__patch_analytical_solution_higher_degree():
    tol=1e-2
    length=48
    height=12
    width=1
    option=3
    nx=27
    ny=13
    nz=3
    load=-20

    Youngs_modulus=100000
    poission_ratio=0.3

    XI_DEGREE = random.randint(2,6)
    ETA_DEGREE = random.randint(2,6)
    NETA_DEGREE = 1

    N = nx
    P = ny
    Q = nz

    C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

    CONTROL_POINTS = C.crtpts_coordinates()
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

    bc_disp = False
    abc_disp = False
    BC = BC_Switcher(CONTROL_POINTS, length, height, width, bc_disp)
    fixed_dof, load_dof, fixed_pts, load_pts = BC.indirect(option)

    reduced_k=np.delete(np.delete(K_G, fixed_dof, 0),fixed_dof , 1)
    reduced_F = apply_BC(F_E, fixed_dof, load_dof, load, abc_disp)


    U = np.linalg.solve(reduced_k, reduced_F)

    for j in fixed_dof:
        U = np.insert(U, j, 0)

    F_E[load_dof] = load

    U_disp=np.concatenate((CONTROL_POINTS[:,:-2], U.reshape(len(CONTROL_POINTS), 3)),axis=1)
    analy_DISP_Y=[]
    Disp_y=[]
    analytical_solution = np.array(exact_displacements(load, length, width, height, Youngs_modulus, poission_ratio, nx, ny,nz))
    #deflection at  point load
    for i in range(len(analytical_solution[:,0])):
        if analytical_solution[i,2]==width/2 and analytical_solution[i,0]==length and analytical_solution[i,1]==height:
            analy_DISP_Y.append(analytical_solution[i,4])
        if U_disp[i,2]==width/2 and U_disp[i,0]==length and U_disp[i,1]==height:   
            Disp_y.append(U_disp[i,4])
    #deflection at  center    
    for i in range(len(analytical_solution[:,0])):
        if analytical_solution[i,2]==0 and analytical_solution[i,0]==length/2 and analytical_solution[i,1]==height/2:
            analy_DISP_Y.append(analytical_solution[i,4])
        if U_disp[i,2]==0 and U_disp[i,0]==length/2 and U_disp[i,1]==height/2:   
            Disp_y.append(U_disp[i,4])
    #deflection at bottom end
    for i in range(len(analytical_solution[:,0])):
        if analytical_solution[i,2]==0 and analytical_solution[i,0]==length and analytical_solution[i,1]==0:
            analy_DISP_Y.append(analytical_solution[i,4])
        if U_disp[i,2]==0 and U_disp[i,0]==length and U_disp[i,1]==0:   
            Disp_y.append(U_disp[i,4])
    expected_output=np.array(analy_DISP_Y)
    output=np.array(Disp_y)
    er=np.all(abs(expected_output-output)<tol)
    o=0
    if er:
        o=1
    assert (o==1) is True
