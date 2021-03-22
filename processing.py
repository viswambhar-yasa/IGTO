#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#---------------------------------------------------------------------------------------#
#A python file where Iso-geometric analysis is performed along with Toplogy optimization
# --------------------------------------------------------------------------------------# 

from inputs import *
from geometry import knot_connectivity, controlpointassembly
from element_routine import assemble, element_routine, apply_BC,Compliance_matrix ,stress_strain_element_routine,gauss_quadrature
from boundary_conditions import BC_Switcher
from optimization import Knearestneighbours, optimality_criteria, Moving_asymptoes
import matplotlib.pyplot as plt
from visulaization import element_density_vis,element_density_slided,element_density_slided1,deformation_plot,mesh_vis
from gridtoVTK import VTK
import time 

#initialization of colour code to plot text on terminal
TYELLOW =  '\033[33;1m' 
TGREEN = '\033[32;1m'
TRED='\033[31;1m'
TBLUE = '\033[34;1m'
ENDC = '\033[m'

#maximum values of penality and optimiality power
pmax = 3.5
gmax = 1

#Degree of knot along xi,eta,neta
XI_DEGREE = 1
ETA_DEGREE = 1
NETA_DEGREE = 1

N = nx
P = ny
Q = nz

# Inputs parameters required for IGA like knot vectors, control points, knot span are generated.
C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

CONTROL_POINTS = C.crtpts_coordinates()

WEIGHTS = CONTROL_POINTS[:, -1]

XI_KNOTVECTOR = C.xi_knotvector()
ETA_KNOTVECTOR = C.eta_knotvector()
NETA_KNOTVECTOR = C.neta_knotvector()
XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()

print('\n')

ncp = N * P * Q  #No of control points 
dof = 3             # Degree of freedom at each point
dofcp = ncp * dof   # Total number of degree of freedom of the structure 
nel = nU * nV * nW  # No of elements
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

#intialization of dimension for Global stiffness matrix, external force vector  and displacements
K_G = np.zeros((dofcp, dofcp))
F_E = np.zeros(dofcp)
U = np.zeros(dofcp)
print('#' * width1)
fmt = '{:^' + str(width1) + '}'
print(TBLUE+fmt.format('Program has started \n')+ENDC)

K_disp = True
#Generation of control point assembly and knot index of each element 
element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                        ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)
span_index = knot_connectivity(N, P, Q, XI_KNOTCONNECTIVITY, ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)

print('$' * width1)
fmt = '{:^' + str(width1) + '}'
IGA_start = time.time() 
KG_ex_time=0
print(TYELLOW+fmt.format('Finite Element Analysis based on ISO-Geometric analysis(NURBS)\n')+ENDC)
#looped over number of element to build global stiffness matrix
for i in range(0, nel):
    KG_start=time.time()  # start of loop 
    el_in = element_indicies[i, :] # element indices obtained from control point assembly contianing the nodes present in the element
    sp_in = span_index[i, :]       # No of knots present with element 
    #Co-ordinates and weights of the control points
    X = CONTROL_POINTS[el_in, 0]   
    Y = CONTROL_POINTS[el_in, 1]
    Z = CONTROL_POINTS[el_in, 2]
    weights = CONTROL_POINTS[el_in, 3]
    #length of the knot in respective direction(xi,eta,neta).
    Uspan = XI_SPAN[sp_in[0], :]
    Vspan = ETA_SPAN[sp_in[1], :]
    Wspan = NETA_SPAN[sp_in[2], :]
    #obtaining element stiffness matrix from element routine.
    K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan, XI_DEGREE,
                                    XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE, NETA_KNOTVECTOR)
    #Assembly of global stiffness matrix from element stiffness matrix
    K_G = assemble(K_G, K_E, el_in, ncp, K_disp)
    K_disp = False
    KG_stop=time.time() # end of loop
    KG_ex_time+=(KG_stop-KG_start) # time taken to build the global stiffness matrix is calculated
print('     Execution time :',KG_ex_time,'\n')
bc_disp = False
abc_disp = True
BC_start=time.time()
#Boundary conditon are obtained from python switch class based on the BC_option
BC = BC_Switcher(CONTROL_POINTS, length, height, width, bc_disp)
fixed_dof, load_dof, fixed_pts, load_pts = BC.indirect(option) #fixed and load indices are obtained

# Boundary conditions are applied by deleting the fixed rows and columns of global stiffness matrix and external force vector.
reduced_k=np.delete(np.delete(K_G, fixed_dof, 0),fixed_dof , 1)
reduced_F = apply_BC(F_E, fixed_dof, load_dof, load,option,abc_disp)

BC_stop=time.time()
BC_ex_time=(BC_stop-BC_start)
print('     Execution time :',BC_ex_time,'\n')

DIS_start=time.time()

# The displacement are calculated from reduced global stiffness matrix and reduced external force vector.
U = np.linalg.solve(reduced_k, reduced_F) # based on equ. 4.15

#Mapping of the displacement along with fixed nodes.
print('Calculating Displacements \n')
for j in fixed_dof:
    U = np.insert(U, j, 0)

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

IGA_stop = time.time() 
IGA_ex_time=(IGA_stop-IGA_start) # time taken to run the IGFEM is calculated
print('Execution time for IGA analysis at 100% volume :',IGA_ex_time,'\n')
energy_stored = np.dot(0.5, U @ F_E)  #strain energy is calculated
print(TRED+'\nThe structure is not optimised and has 100% volume \n'+ENDC)
print('The strain energy of the structure               :', energy_stored)
print('$' * width1)
CP = CONTROL_POINTS[:, :-2]
if mesh_disp:
    #plotting the deformed and undeformed 3D structure. 
    mesh_vis(CP, New_control_points, nx, ny, nz,optimizer)



# To implement modified SIMP method equ. 4.21, we require Young's modulus of void and solid
Emin = 1e-09  # Young's modulus of void (intializing a value almost equal to 0 to avoid null point errors)
E0 = 1 # Young's mosulus of solid
#intializing the dimensions of global stiffness matrix and external force vector for toplogy optimization
K_G = np.zeros((dofcp, dofcp))
F_E = np.zeros(dofcp)
U = np.zeros(dofcp)
CC = []
ii = []
VV = []
oc_disp = True
bc_disp = True
fil_disp = True
max_iterations = 250
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
#initialzing dimension and initial values of variables used in toplogy optimization
density_basis = np.zeros(nel)
ele_den_filter = np.zeros(nel)
element_density = np.ones(nel) * volume_frac # density of each element (initially taken as one)
density_basis_dx = np.zeros(nel)
dcompliance = np.zeros(nel)  # change in compliance w.r.t element density (initially taken as zero)
compliance = 0 #intial compliance of the structure
# performing sensitivity analysis i.e giving weights to the respective elements.
nfilter = int(nel * ((2 * (np.ceil(rmin) - 1) + 1) ** 2))
filter_N = np.zeros(((XI_DEGREE + 1) * (ETA_DEGREE + 1) * (NETA_DEGREE + 1), nel))
#Weight factor are obtained from the below function
H, DH = Knearestneighbours(rmin, nU, nV, nW) 
loop = 0
change = 1
g = 1

if penal==0:
    penal=max(15*((1-poission_ratio)/(7-5*poission_ratio)),(3/2)*((1-poission_ratio)/(1-2*poission_ratio)))

#intializing dimension and initial value of the variables.
Xmin = np.zeros(nel) # Minimum of the element density i.e 0
Xmax = np.ones(nel)  # Maximum of the element density i.e 1
#Lower and Upper values are used in MMA method
Lower = Xmin   
Upper = Xmin
E1 = element_density # values of element density in previous iteration k-1
E2 = element_density # value of element density in previous iteration k-2
Total_time=0
OC_residual=0
# loop run until the termination condition is satisified i.e change in compliance from previous iteration 
while change > 0.01:
    #intializing the dimensions of global stiffness matrix, external force vector, compliance,change in compliance for toplogy optimization in each iteration.
    K_G = np.zeros((dofcp, dofcp))
    F_E = np.zeros(dofcp)
    U = np.zeros(dofcp)
    dcompliance = np.zeros(nel)
    compliance = 0
    # nodal density are calculated based on equ. 4.20
    node_density = np.ones((nel, (XI_DEGREE + 1) * (ETA_DEGREE + 1) * (NETA_DEGREE + 1)))
    for h in range(nel):
        node_density[h, :] = node_density[h, :] * element_density[h]  #Based on equ. 4.20
    #looped over number of elements
    for i in range(0, nel):
        el_in = element_indicies[i, :]  # control point assembly
        sp_in = span_index[i, :]        #span each knot along xi,eta,neta
        #co-ordinates of the control points
        X = CONTROL_POINTS[el_in, 0]
        Y = CONTROL_POINTS[el_in, 1]
        Z = CONTROL_POINTS[el_in, 2]
        weights = CONTROL_POINTS[el_in, 3]
        Uspan = XI_SPAN[sp_in[0], :]
        Vspan = ETA_SPAN[sp_in[1], :]
        Wspan = NETA_SPAN[sp_in[2], :]
        #Obtaining element stiffness matrix
        K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan,
                                        XI_DEGREE, XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE,
                                        NETA_KNOTVECTOR)
        #calculating element density based on equ. 4.20                                
        element_density[i] = np.dot(node_density[i, :], NURBS)
        filter_N[:, i] = NURBS
        #Implementation of Modified SIMP method
        density_basis[i] = Emin + (element_density[i] ** penal) * (E0 - Emin) #Based on equ. 4.21
        # Global stiffness matrix is calculated 
        K_E = density_basis[i] * K_E  
        K_G = assemble(K_G, K_E, el_in, ncp) #based on equ. 4.17
    
    #Boundary condition are obtained 
    BC = BC_Switcher(CONTROL_POINTS, length, height, width, bc_disp)
    fixed_dof, load_dof, fixed_pts, load_pts = BC.indirect(option)
    #Boundary conditions are applied 
    reduced_k=np.delete(np.delete(K_G, fixed_dof, 0),fixed_dof , 1)
    reduced_F = apply_BC(F_E, fixed_dof, load_dof, load,option)
    # Displacements ara calculated 
    U = np.matmul(np.linalg.inv(reduced_k), reduced_F) #based on equ. 4.17
    #Displacements are mapped
    for j in fixed_dof:
        U = np.insert(U, j, 0)

    if option==4:
        F_E[load_dof[0]] = load
        F_E[load_dof[1]] = -load
    else:
        F_E[load_dof] = load
    # Untill now, The displacements of the optimized structure are calculated whose dispalcements are required to compute compliance and change in 
    # compliance w.r.t element density. 

    #looped over number of no of element to calculate compliance and dcompliance 
    for k in range(0, nel):
        loop_start=time.time() #start of the iteration
        el_in = element_indicies[k, :] #control point assembly
        sp_in = span_index[k, :]      #knot span along xi,eta and neta directon
        #co-ordinates of the control points
        X = CONTROL_POINTS[el_in, 0]
        Y = CONTROL_POINTS[el_in, 1]
        Z = CONTROL_POINTS[el_in, 2]
        weights = CONTROL_POINTS[el_in, 3]
        Uspan = XI_SPAN[sp_in[0], :]
        Vspan = ETA_SPAN[sp_in[1], :]
        Wspan = NETA_SPAN[sp_in[2], :]
        #obtained element stiffness matrix
        K_E, NURBS, R,B = element_routine(X, Y, Z, weights, Youngs_modulus, poission_ratio, Uspan, Vspan, Wspan,
                                        XI_DEGREE, XI_KNOTVECTOR, ETA_DEGREE, ETA_KNOTVECTOR, NETA_DEGREE,
                                        NETA_KNOTVECTOR)
        #Calculating element density based on equ. 4.20                                
        element_density[k] = np.dot(node_density[k, :], NURBS)
        # Implemenation of SIMP method based on equ. 4.21
        density_basis[k] = Emin + (element_density[k] ** penal) * (E0 - Emin)
        #derivative of equ. 4.21 w.r.t to element density
        density_basis_dx[k] = -penal * (element_density[k] ** (penal - 1)) * (E0 - Emin)
        dof_index = np.sort(np.concatenate((el_in * dof, dof * el_in + 1, dof * el_in + 2)))
        #Calculating compliance and change in compliance
        compliance += np.transpose(U[dof_index]) @ (density_basis[k] * K_E) @ U[dof_index]        #Based on equ. 4.23 Objective function
        dcompliance[k] = np.transpose(U[dof_index]) @ (density_basis_dx[k] * K_E) @ U[dof_index]  #Based on equ. 4.30
    #Initializing the change in volume fraction V0/V
    dv = np.ones(nel)

    element_density = np.round(element_density, 5)
    old_el_density = element_density # previous element density are stored 

    constrain = nel * volume_frac
    #To solve constrained based probelm, we used methods like penality method(OC) and active set strategy(MMA) to satisfy KKT conditions
    if optimizer == 'OC': #Optimality criteria
        
        #dcompliance = (1 / np.maximum(1e-3, element_density) * DH) * np.matmul(H, (element_density * dcompliance))
        #implementation of sensitivity analysis
        dcompliance=np.matmul(H,(dcompliance/DH))  #Based on equ. 4.34
        dv = np.matmul(H, dv / DH)                  #based on equ. 4.35
        
        element_density_updated,OC_residual = optimality_criteria(dcompliance, dv, element_density, constrain,OC_residual, H, DH, oc_disp,
                                                      g)
        element_density = element_density_updated
    if optimizer == 'MMA': #Method of moving asymptoes
        #Sensitivity analysis 
        dcompliance = np.matmul(H, (dcompliance / DH))  #Based on 4.34
        dv = np.matmul(H, dv / DH)                      #based on 4.35

        v = (np.sum(element_density) / nel) - volume_frac

        dv_dx = (dv / (volume_frac * nel))

        dfun0 = dcompliance
        dcon = dv_dx
        f0 = compliance
        c0 = v
        element_density_updated, Lower, Upper = Moving_asymptoes(dfun0, dcon, f0, c0, element_density, E1, E2, Lower,
                                                                 Upper, loop, nel, 1, Xmin, Xmax, True)
        
        E2 = E1.copy()
        E1 = element_density.copy()
        element_density = np.round(element_density_updated[:, 0], 4)
    #Compliance and volume after each iteration are stored for analysis
    CC.append(compliance)
    ii.append(loop)
    VV.append(np.mean(element_density))
    #Change in the element density is calculated
    change = np.linalg.norm(element_density - old_el_density)

    if loop == 0:
        initial_compliance = compliance
    loop += 1
    print(TYELLOW+'     Iteration: %d       p :%f      optimizer:%s      OBJ=%f6       Vol=%f       CHG=%f6      ' % (
        loop, penal, optimizer, compliance, np.mean(element_density), change)+ENDC)
    print('#' * width1)
    loop_stop=time.time()
    loop_ex_time=loop_stop-loop_start
    Total_time+=loop_ex_time #time taken to run the toplogy optimization loop is calculated
     
    print('       Execution time of this iteration :',loop_ex_time,'sec          Total execution Time:',Total_time,'sec') 
    print('=' * width1)
    print('\n')
    CPS = CONTROL_POINTS[:, :-2]
    #Plotting of the optimized structure using PYVISTA library
    if (change<0.01 or (loop in [1,10,20,40,60,80,100,120,140,160,180,200,220,240,250])) and iterative_display is True:
        ns=nU*10
        element_density_slided1(loop,CP,nx,ny,nz,element_density,optimizer)
        #deformation_plot(loop,CP,U,nx,ny,nz,element_density,optimizer)
    if loop >= max_iterations: #termination condition 
        print('~' * width1)
        fmt = '{:^' + str(width1) + '}'
        print(fmt.format('Maximum iteration has been reached'))
        print('~' * width1)
        break
print('+' * width1)
#Measure of discretness - Calculated to show the efficiency of the method's abitlity to penalize the element desnsity as 0 or 1
Mnd = np.dot(4 * element_density, (1 - element_density)) / nel * 100
print('                                        Measure of discreteness: ', Mnd)

element_density = np.round(element_density, 1) # FINAL OPTIMIZED STRUCTURE IS OBTAINED 

print('\n                               Final Volume of the optimised structure  :', np.mean(element_density))
print('\n                               Mass of the beam at 100% volume          :',
      (length * width * height) * density)
print('\n                               Optimized mass of the beam               :',
      (length * width * height) * density * np.mean(element_density))
energy_stored = np.dot(0.5, np.matmul(U.transpose(), np.matmul(U, K_G))) # Strain energy of the optimised structure is calculated.
print('\n                               Initial Compliance without optimization  :', initial_compliance)
final_compliance = np.dot(U, F_E)
print('\n                               Final Compliance with optimization       :', final_compliance, '\n')

#VTK file is genetated
VTK(CONTROL_POINTS,element_density,nU,nV,nW,"Cantilever_beam")


#element_density_vis(loop, CPS, nx, ny, nz, element_density,optimizer)
#element_density_slided(i,CP,nx,ny,nz,element_density,optimizer,ns,'z','xz',False)

