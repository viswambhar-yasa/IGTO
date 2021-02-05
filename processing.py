from Preprocessing import *
from geometry import knot_connectivity, controlpointassembly
from element_routine import assemble, element_routine, apply_BC,Compliance_matrix ,stress_strain_element_routine,gauss_quadrature
from boundary_conditions import BC_Switcher
from optimization import Knearestneighbours, optimality_criteria, Heavyside_filter, Moving_asymptoes
import matplotlib.pyplot as plt
from visulaization import element_density_vis,element_density_slided,element_density_slided1,deformation_plot,mesh_vis
from gridtoVTK import VTK
import time 

'''
length=48
height=12
width=1
option=2
nx=2 
ny=3
nz=3
density=7850
volume_frac=0.5
pmax=5
gmax=1
rmin=1.5
load=-20
'''
TYELLOW =  '\033[33;1m' 
TGREEN = '\033[32;1m'
TRED='\033[31;1m'
TBLUE = '\033[34;1m'
ENDC = '\033[m'

pmax = 3.5
gmax = 1

XI_DEGREE = 1
ETA_DEGREE = 1
NETA_DEGREE = 1

N = nx
P = ny
Q = nz

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

print('\n')

ncp = N * P * Q
dof = 3
dofcp = ncp * dof
nel = nU * nV * nW
width1 = 120
print('#' * width1)
fmt = '{:^' + str(width1) + '}'
print(TBLUE+fmt.format('Dimensions of the structure \n')+ENDC)
print('                   Length       :', length, '               Height     :', height,
      '               Width        :', width, '\n')

print('                   Xi degree    :', XI_DEGREE, '                Eta degree :', ETA_DEGREE,
      '                Neta degree  :', NETA_DEGREE, '\n')
print('                   NX           :', N - XI_DEGREE, '                NY         :', P - ETA_DEGREE,
      '                NZ           :', Q - NETA_DEGREE, '\n')
print(TGREEN+'Number of degrees of freedom :', dofcp, '\n')

print('Number of Elements:', nel, '\n')

print('No of control points in each element:', (XI_DEGREE + 1) * (ETA_DEGREE + 1) * (NETA_DEGREE + 1), '\n'+ENDC)

print('>' * width1)
print('Length of the knot vector in respective direction \n')
print('XI Vector        :', list(XI_KNOTVECTOR), '\nETA vector       :', list(ETA_KNOTVECTOR), '\nNETA vector      :',
      list(NETA_KNOTVECTOR), '\n')
print('<' * width1)
K_G = np.zeros((dofcp, dofcp))
F_E = np.zeros(dofcp)
U = np.zeros(dofcp)
print('#' * width1)
fmt = '{:^' + str(width1) + '}'
print(TBLUE+fmt.format('Program has started \n')+ENDC)

K_disp = True
element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                        ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
span_index = knot_connectivity(N, P, Q, XI_KNOTCONNECTIVITY, ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
print('$' * width1)
fmt = '{:^' + str(width1) + '}'
IGA_start = time.time() 
KG_ex_time=0
print(TYELLOW+fmt.format('Finite Element Analysis based on ISO-Geometric analysis(NURBS)\n')+ENDC)
for i in range(0, nel):
    KG_start=time.time()  
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
    KG_stop=time.time()
    KG_ex_time+=(KG_stop-KG_start)
print('     Execution time :',KG_ex_time,'\n')
bc_disp = False
abc_disp = True
BC_start=time.time()
BC = BC_Switcher(CONTROL_POINTS, length, height, width, bc_disp)
fixed_dof, load_dof, fixed_pts, load_pts = BC.indirect(option)

reduced_k=np.delete(np.delete(K_G, fixed_dof, 0),fixed_dof , 1)
reduced_F = apply_BC(F_E, fixed_dof, load_dof, load,option,abc_disp)

BC_stop=time.time()
BC_ex_time=(BC_stop-BC_start)
print('     Execution time :',BC_ex_time,'\n')

DIS_start=time.time()

U = np.linalg.solve(reduced_k, reduced_F)
#print(U)
print('Calculating Displacements \n')
for j in fixed_dof:
    U = np.insert(U, j, 0)
#print(U)
if option==4:
    F_E[load_dof[0]] = load
    F_E[load_dof[1]] = -load
else:
    F_E[load_dof] = load

DIS_stop=time.time()
DIS_ex_time=(DIS_stop-DIS_start)
print('     Execution time :',DIS_ex_time,'\n')
U_new = np.array((U.reshape(len(CONTROL_POINTS), 3)), dtype='float64')
print('Mapping Displacements \n')
New_control_points = CONTROL_POINTS[:,:-2]+ U.reshape(len(CONTROL_POINTS), 3)
UX = U_new[:, 0]
UY = U_new[:, 1]
UZ = U_new[:, 2]
#print(np.concatenate((CONTROL_POINTS[:,:-2], U.reshape(len(CONTROL_POINTS), 3)),axis=1))
#print(U_new)
#print(New_control_points)

gauss_points,gauss_weights=gauss_quadrature(XI_DEGREE,ETA_DEGREE,NETA_DEGREE)
gauss_points_strains=np.zeros((nel,len(gauss_points),6))
gauss_point_stresses=np.zeros((nel,len(gauss_points),6))
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
    gauss_points_strains[i,:,:],gauss_point_stresses[i,:,:],gauss_point_disp[i,:,:] = stress_strain_element_routine(X, Y, Z, weights,U_el, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                    XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)                              
    #print(gauss_point_stresses[i,:,0])
strain_xx=strain_yy=strain_zz=strain_xy=strain_yz=strain_xz=np.zeros(len(CONTROL_POINTS[:,0]))
sigma_xx=sigma_yy=sigma_zz=sigma_xy=sigma_yz=sigma_xz=np.zeros(len(CONTROL_POINTS[:,0]))
Disp_x=Disp_y=Disp_z=np.zeros(len(CONTROL_POINTS[:,0]))
print(gauss_point_disp)
for i in range(0,nel):
    el_in = element_indicies[i, :]

    sigma_xx[el_in]=gauss_point_stresses[i,:,0]
    #sigma_yy[el_in]=gauss_point_stresses[i,:,1]
    #sigma_zz[el_in]=gauss_point_stresses[i,:,2]
    #sigma_xy[el_in]=gauss_point_stresses[i,:,3]
    #sigma_yz[el_in]=gauss_point_stresses[i,:,4]
    #sigma_xz[el_in]=gauss_point_stresses[i,:,5]
    #strain_xx[el_in]=gauss_points_strains[i,:,0]
    #strain_yy[el_in]=gauss_points_strains[i,:,1]
    #strain_zz[el_in]=gauss_points_strains[i,:,2]
    #strain_xy[el_in]=gauss_points_strains[i,:,3]
    #strain_yz[el_in]=gauss_points_strains[i,:,4]
    #strain_xz[el_in]=gauss_points_strains[i,:,5]
    Disp_x[el_in]=gauss_point_disp[i,:,0]
    Disp_y[el_in]=gauss_point_disp[i,:,1]
    Disp_z[el_in]=gauss_point_disp[i,:,2]
#print('gauss_point_disp' ,gauss_point_disp[0,:,0])
 
IGA_stop = time.time() 
IGA_ex_time=(IGA_stop-IGA_start)
print('Execution time for IGA analysis at 100% volume :',IGA_ex_time,'\n')
energy_stored = np.dot(0.5, U @ F_E)
print(TRED+'\nThe structure is not optimised and has 100% volume \n'+ENDC)
print('The strain energy of the structure               :', energy_stored)
print('$' * width1)
CP = CONTROL_POINTS[:, :-2]
'''
mesh_vis(CP, New_control_points, nx, ny, nz,optimizer)
'''

Emin = 1e-09
E0 = 1
K_G = np.zeros((dofcp, dofcp))
F_E = np.zeros(dofcp)
U = np.zeros(dofcp)
print('+' * width1)
fmt = '{:^' + str(width1) + '}'
print(TYELLOW+fmt.format('Structural Topology optimization using IGA \n')+ENDC)
print(TBLUE+'\n                                            Optimization has started \n'+ENDC)

print('                                     Density of the material       :', density)
print('                                     Youngs Modulus                :', Youngs_modulus)
print('                                     Poission ratio                :', poission_ratio, '\n')
fmt = '{:^' + str(width1) + '}'
print(TRED+fmt.format('The percentage of volume which has to remain after optimization \n'))
fmt = '{:^' + str(width1) + '}'
print(fmt.format( volume_frac)+ENDC)
print('+' * width1)
density_basis = np.zeros(nel)
ele_den_filter = np.zeros(nel)
element_density = np.ones(nel) * volume_frac
density_basis_dx = np.zeros(nel)
dcompliance = np.zeros(nel)
compliance = 0
nfilter = int(nel * ((2 * (np.ceil(rmin) - 1) + 1) ** 2))
H, DH = Knearestneighbours(rmin, nU, nV, nW)
loop = 0
change = 1
g = 1
CC = []
ii = []
VV = []
# penal=max(15*((1-poission_ratio)/(7-5*poission_ratio)),(3/2)*((1-poission_ratio)/(1-2*poission_ratio)))
oc_disp = True
bc_disp = True
fil_disp = True
max_iterations = 250
beta =1
filter_N = np.zeros(((XI_DEGREE + 1) * (ETA_DEGREE + 1) * (NETA_DEGREE + 1), nel))
Xmin = np.zeros(nel)
Xmax = np.ones(nel)
Lower = Xmin
Upper = Xmin
E1 = element_density
E2 = element_density
Total_time=0
while change > 0.01:
    K_G = np.zeros((dofcp, dofcp))
    F_E = np.zeros(dofcp)
    U = np.zeros(dofcp)
    dcompliance = np.zeros(nel)
    compliance = 0
    node_density = np.ones((nel, (XI_DEGREE + 1) * (ETA_DEGREE + 1) * (NETA_DEGREE + 1)))
    for h in range(nel):
        node_density[h, :] = node_density[h, :] * element_density[h]
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

        K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan,
                                        XI_DEGREE, XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE,
                                        NETA_KNOTVECTOR)
        element_density[i] = np.dot(node_density[i, :], NURBS)
        filter_N[:, i] = NURBS
        density_basis[i] = Emin + (element_density[i] ** penal) * (E0 - Emin)
        K_E = density_basis[i] * K_E
        K_G = assemble(K_G, K_E, el_in, ncp)

    BC = BC_Switcher(CONTROL_POINTS, length, height, width, bc_disp)
    fixed_dof, load_dof, fixed_pts, load_pts = BC.indirect(option)
    
    reduced_k=np.delete(np.delete(K_G, fixed_dof, 0),fixed_dof , 1)
    reduced_F = apply_BC(F_E, fixed_dof, load_dof, load,option)
    U = np.matmul(np.linalg.inv(reduced_k), reduced_F)

    for j in fixed_dof:
        U = np.insert(U, j, 0)

    if option==4:
        F_E[load_dof[0]] = load
        F_E[load_dof[1]] = -load
    else:
        F_E[load_dof] = load
    # print(U@F_E)

    for k in range(0, nel):
        loop_start=time.time()
        el_in = element_indicies[k, :]
        sp_in = span_index[k, :]
        X = CONTROL_POINTS[el_in, 0]
        Y = CONTROL_POINTS[el_in, 1]
        Z = CONTROL_POINTS[el_in, 2]
        weights = CONTROL_POINTS[el_in, 3]
        Uspan = XI_SPAN[sp_in[0], :]
        Vspan = ETA_SPAN[sp_in[1], :]
        Wspan = NETA_SPAN[sp_in[2], :]

        K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan,
                                        XI_DEGREE, XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE,
                                        NETA_KNOTVECTOR)
        element_density[k] = np.dot(node_density[k, :], NURBS)

        density_basis[k] = Emin + (element_density[k] ** penal) * (E0 - Emin)
        density_basis_dx[k] = -penal * (element_density[k] ** (penal - 1)) * (E0 - Emin)
        dof_index = np.sort(np.concatenate((el_in * dof, dof * el_in + 1, dof * el_in + 2)))
        compliance += np.transpose(U[dof_index]) @ (density_basis[k] * K_E) @ U[dof_index]
        dcompliance[k] = np.transpose(U[dof_index]) @ (density_basis_dx[k] * K_E) @ U[dof_index]
    # compliance+=np.transpose(U)@(K_G)@U
    dv = np.ones(nel)
    # dcompliance=(1/np.maximum(1e-3,element_density)*DH)*np.matmul(H,(element_density*dcompliance))
    # delement_density=beta*np.exp(-beta*element_density)+np.exp(-beta)
    # dcompliance=np.matmul(H,(delement_density*dcompliance/DH))
    # dcompliance=np.matmul(H,(dcompliance/DH))
    # dv=np.matmul(H,dv/DH)
    # print(dcompliance)
    element_density = np.round(element_density, 5)
    old_el_density = element_density
    # np.append(X,element_density)
    # element_density=np.matmul(H,element_density/DH)
    constrain = nel * volume_frac
    if optimizer == 'OC':
        
        dcompliance = (1 / np.maximum(1e-3, element_density) * DH) * np.matmul(H, (element_density * dcompliance))
        dv = np.matmul(H, dv / DH)
        
        element_density_updated = optimality_criteria(dcompliance, dv, element_density, constrain, H, DH, beta, oc_disp,
                                                      g)
        element_density = element_density_updated
    if optimizer == 'MMA':
        dcompliance = np.matmul(H, (dcompliance / DH))
        dv = np.matmul(H, dv / DH)

        v = (np.sum(element_density) / nel) - volume_frac
        # print(v)

        dv_dx = (dv / (volume_frac * nel))
        # print(dcompliance)
        # print(dv)
        dfun0 = dcompliance
        dcon = dv_dx
        f0 = compliance
        c0 = v
        element_density_updated, Lower, Upper = Moving_asymptoes(dfun0, dcon, f0, c0, element_density, E1, E2, Lower,
                                                                 Upper, loop, nel, 1, Xmin, Xmax, True)
        E2 = E1.copy()
        E1 = element_density.copy()
        element_density = np.round(element_density_updated[:, 0], 4)
    # Heavyside filter
    # filter_density=Heavyside_filter(element_density,beta,fil_disp)
    # print(element_density)
    # element_density=filter_density
    if (loop % 50) == 0:
        beta *= 2
    CC.append(compliance)
    ii.append(loop)
    VV.append(np.mean(element_density))
    change = np.linalg.norm(element_density - old_el_density)
    # if optimizer=='MMA':
    #    element_density=np.matmul(H,element_density/DH)
    if loop == 0:
        initial_compliance = compliance
    loop += 1
    print(TYELLOW+'     Iteration: %d       p :%f      optimizer:%s      OBJ=%f6       Vol=%f       CHG=%f6      ' % (
        loop, penal, optimizer, compliance, np.mean(element_density), change)+ENDC)
    print('#' * width1)
    loop_stop=time.time()
    loop_ex_time=loop_stop-loop_start
    Total_time+=loop_ex_time
     
    print('       Execution time of this iteration :',loop_ex_time,'sec          Total execution Time:',Total_time,'sec') 
    print('=' * width1)
    print('\n')
    CPS = CONTROL_POINTS[:, :-2]
    if iterative_display or change<0.01 or ((loop%10)==0):
        ns=nU*10
        deformation_plot(loop,CP,U,nx,ny,nz,element_density,optimizer)
        #element_density_slided(i,CP,nx,ny,nz,element_density,optimizer,ns,'z','xz')
        element_density_slided1(loop,CP,nx,ny,nz,element_density,optimizer)
    if loop >= max_iterations:
        print('~' * width1)
        fmt = '{:^' + str(width1) + '}'
        print(fmt.format('Maximum iteration has been reached'))
        print('~' * width1)
        break
print('+' * width1)

Mnd = np.dot(4 * element_density, (1 - element_density)) / nel * 100
print('                                        Measure of discreteness: ', Mnd)

element_density = np.round(element_density, 1)
# element_density[element_density<0.1] = 0
# element_density[element_density>=volume_frac] = 1

print('\n                               Final Volume of the optimised structure  :', np.mean(element_density))
print('\n                               Mass of the beam at 100% volume          :',
      (length * width * height) * density)
print('\n                               Optimized mass of the beam               :',
      (length * width * height) * density * np.mean(element_density))
energy_stored = np.dot(0.5, np.matmul(U.transpose(), np.matmul(U, K_G)))
print('\n                               Initial Compliance without optimization  :', initial_compliance)
final_compliance = np.dot(U, F_E)
print('\n                               Final Compliance with optimization       :', final_compliance, '\n')


VTK(CONTROL_POINTS,element_density,nU,nV,nW,"Cantilever_beam")
#element_density_vis(loop, CPS, nx, ny, nz, element_density,optimizer)
#element_density_slided(i,CP,nx,ny,nz,element_density,optimizer,ns,'z','xz',False)

