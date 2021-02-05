from test_preprocessing import *
from geometry import knot_connectivity, controlpointassembly
from element_routine import assemble, element_routine, apply_BC,Compliance_matrix ,stress_strain_element_routine,gauss_quadrature
from boundary_conditions import BC_Switcher
from analytical_solution import exact_displacements
import math as m
import time 
import pytest
import random


def test__patch_constant_strain():
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
    #print(CONTROL_POINTS)
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

    F_residual=np.matmul(K_G,U)
    disp_dof=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26])
    disp_index=np.sort(np.concatenate((disp_dof*3,disp_dof*3+1,disp_dof*3+2)))
    load_dof=np.array([13])
    #load_index=np.array([13*3,13*3+1,13*3+2])
    load_index=np.sort(np.concatenate((load_dof*3,load_dof*3+1,load_dof*3+2)))
    U[load_index]=(F_E[load_index]-F_residual[load_index])/(K_G[load_index,load_index])
    #print(U)
    gauss_points,gauss_weights=gauss_quadrature(XI_DEGREE,ETA_DEGREE,NETA_DEGREE)
    gauss_points_strains=np.zeros((nel,len(gauss_points),6))
    gauss_point_stresses=np.zeros((nel,len(gauss_points),6))
    nodal_stress=np.zeros((nel,8,6))
    nodal_strain=np.zeros((nel,8,6))
    gauss_point_disp=np.zeros((nel,len(gauss_points),3))
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
        U_el=U[u_index].T
        gauss_points_strains[i,:,:],gauss_point_stresses[i,:,:],gauss_point_disp[i,:,:],nodal_stress[i,:,:],nodal_strain[i,:,:] = stress_strain_element_routine(X, Y, Z, weights,U_el, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                        XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)                              
    print(np.round(gauss_points_strains,3))
    assert gauss_points_strains==0.001 or  gauss_points_strains==0


def test__comparing_with_analytical_solution():
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
    if er==True:
        o=1
    assert (o==1) is True

def test__comparing_with_analytical_solution_higher_degree():
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
    if er==True:
        o=1
    assert (o==1) is True

''' 
disp_index=np.array([21,22,23])
U[disp_index]=np.matmul(np.array([1,1,1]),rotation_matrix)
WEIGHTS = CONTROL_POINTS[:, -1]
print(U)
residual_F=np.matmul(K_G,U)
index=np.array([0,1,2,3,4,5,6,7])
fixed_index=np.sort(np.concatenate((np.array([0,1,12,13]),index*3+2)))
F=F_E-residual_F

reduced_k=np.delete(np.delete(K_G, fixed_index, 0),fixed_index , 1)
reduced_F = np.delete(F,fixed_index)
U = np.linalg.solve(reduced_k, reduced_F)
for j in fixed_index:
    U = np.insert(U, j, 0)
print(np.round(New_control_points,5))
print(ctrl_pts+U.reshape(len(ctrl_pts),3))*

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
print(CONTROL_POINTS)
CONTROL_POINTS[13,:]=np.array([ 0.6,0.4,0.3,1.,13.])
#print(CONTROL_POINTS)
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
F_residual=np.matmul(K_G,U)
disp_dof=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26])
disp_index=np.sort(np.concatenate((disp_dof*3,disp_dof*3+1,disp_dof*3+2)))
load_dof=np.array([13])
#load_index=np.array([13*3,13*3+1,13*3+2])
load_index=np.sort(np.concatenate((load_dof*3,load_dof*3+1,load_dof*3+2)))
U[load_index]=(F_E[load_index]-F_residual[load_index])/(K_G[load_index,load_index])
#print(U)
gauss_points,gauss_weights=gauss_quadrature(XI_DEGREE,ETA_DEGREE,NETA_DEGREE)
gauss_points_strains=np.zeros((nel,len(gauss_points),6))
gauss_point_stresses=np.zeros((nel,len(gauss_points),6))
nodal_stress=np.zeros((nel,8,6))
nodal_strain=np.zeros((nel,8,6))
gauss_point_disp=np.zeros((nel,len(gauss_points),3))
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
    U_el=U[u_index].T
    gauss_points_strains[i,:,:],gauss_point_stresses[i,:,:],gauss_point_disp[i,:,:],nodal_stress[i,:,:],nodal_strain[i,:,:] = stress_strain_element_routine(X, Y, Z, weights,U_el, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                    XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)                              
print(np.round(gauss_points_strains,3))

'''
tol=1e-2
length=2
height=0.1
width=2
option=1
nx=21
ny=2
nz=21
load=4
Youngs_modulus=1000000
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
print(CONTROL_POINTS)
print(fixed_pts,fixed_dof)
reduced_k=np.delete(np.delete(K_G, fixed_dof, 0),fixed_dof , 1)
reduced_F = apply_BC(F_E, fixed_dof, load_dof, load, abc_disp)
U = np.linalg.solve(reduced_k, reduced_F)
for j in fixed_dof:
    U = np.insert(U, j, 0)
F_E[load_dof] = load

print(U[load_dof]) 
U_disp=np.concatenate((CONTROL_POINTS[:,:-2], U.reshape(len(CONTROL_POINTS), 3)),axis=1)
analytical_solution=(0.1267*load*(length**2))/(Youngs_modulus*(height**3)) 
print(analytical_solution)