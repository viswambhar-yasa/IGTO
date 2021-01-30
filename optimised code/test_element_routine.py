from Preprocessing import *
from element_routine import gauss_quadrature,unittoparametric,jacobian,assemble,element_routine,apply_BC 
from numpy import array,array_equiv,zeros,linalg,sin,cos,radians,transpose,min,round
from Inputs import *
from geometry import knot_connectivity,controlpointassembly
from boundary_conditions import BC_Switcher
import pytest



def test__gauss_quadrature_1point_true():
    gausspoints,weights=gauss_quadrature(0,0,0)
    assert( array_equiv(gausspoints,array([0.0,0.0,0.0])) and array_equiv(weights,array([8.0])) )is True

def test__unit_to_parametric_space_true():
    '''
    value calculated anaytically 

    '''
    output=float(round(unittoparametric(-0.5773,[0,0.2]),5))
    assert (output==float(0.04227)) is True



def test__Jacobian_patch_testing_rotation():
    '''
    rotation of the jacbian matrix using Bunge convenstion should change the determinat as it is volume conserving
    '''
    
    X=unittoparametric(-0.5773,[0,1])
    E=unittoparametric(0.5773,[0,1])
    N=unittoparametric(-0.5773,[0,1])
    Wspa=Vspa=Uspa=array([0,1])
    control_points=array([[0., 0., 0. ,1.],
                    [1., 0., 0. ,1.],
                    [0., 1., 0., 1.],
                    [1., 1., 0., 1.],
                    [0., 0., 1., 1.],
                    [1., 0., 1., 1.],
                    [0., 1., 1., 1.],
                    [1., 1., 1., 1.]])
    xx=control_points[:,0]
    yy=control_points[:,1]
    zz=control_points[:,2]
    weigt=control_points[:,3]
    xdef=1
    XKV=[0, 0, 1, 1]
    jacobian1,N1,N2,N3,N4=jacobian(X,E,N,Uspa,Vspa,Wspa,xx,yy,zz,weigt,xdef,XKV,xdef,XKV,xdef,XKV)
    output=zeros(360)
    for i in range(0,360):
        a=radians(i)
        Q1=array([[cos(a),sin(a) ,0],[-sin(a),cos(a) ,0],[0,0,1]])
        Q2=array([[cos(a),0,-sin(a)],[0,1,0],[sin(a),0,cos(a)]])
        Q3=array([[cos(a),sin(a) ,0],[-sin(a),cos(a) ,0],[0,0,1]])
        Q=Q1@Q2@Q3
        Q_T=transpose(Q)
        output[i]=round(linalg.det(Q_T@jacobian1@Q),10)
    actual_det_jacobian=round(linalg.det(jacobian1),10)
    assert all(output==actual_det_jacobian) is True  

def test__stiffness_matrix_singularity():
    '''
    Stiffness matrix shold be positive definit

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
    
    eigen_values,eigen_vectors= linalg.eig(K_G)
    eigen_values=round(eigen_values,8)
    assert (all(eigen_values>=0)) is True


def test__rigid_body_motion():
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

    eigen_values,eigen_vectors= linalg.eig(K_G)
    eigen_values=round(eigen_values,8)
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
    eigen_values,eigen_vectors= linalg.eig(reduced_k)
    eigen_values=round(eigen_values,8)
    print(eigen_values)
    constrained_rigid_body_modes=len([eigen_values[eigen_values==0]])
    #Too few zero eigenvalues is an indication of an element lacking the capability of rigid body motion without strain.
    #Translation along y axis Cantilever beam
    assert (unconstrained_rigid_body_modes==6) and (constrained_rigid_body_modes==1) is True

length=10
height=10
width=0.1
nx=3
ny=3
nz=2
XI_DEGREE=1
ETA_DEGREE=1
NETA_DEGREE=1
Youngs_modulus=1000
poission_ratio=0.25
N = nx
P = ny
Q = nz
C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)
CONTROL_POINTS = C.crtpts_coordinates()
U=[]
'''
for i  in range(len(CONTROL_POINTS[:,0])):
    x1=CONTROL_POINTS[i,0]
    y1=CONTROL_POINTS[i,1]
    z1=CONTROL_POINTS[i,2]
    Ux=1e-3*(2*x1+y1+z1)/2
    Uy=1e-3*(x1+2*y1+z1)/2
    Uz=1e-3*(x1+y1+2*z1)/2
    U.append(Ux)
    U.append(Uy)
    U.append(Uz)
'''
#U=np.concatenate(Ux.T,np.concatenate(Uy,Uz))
#U[(13*3)]=U[(13*3+1)]=U[(13*3+2)]=0
#load_index=np.array([13])
#loading_index=np.sort(np.concatenate((load_index*3,3*load_index+1,load_index*3+2)))
#U=np.array(U)
#CONTROL_POINTS[13,:]=np.array([ 0.6,0.4,0.3,1.,13.])

CONTROL_POINTS[1,:]=np.array([ 4,0,0,1.,1.])
CONTROL_POINTS[10,:]=np.array([ 4,0,0,1.,10.])
CONTROL_POINTS[7,:]=np.array([ 4.2,10,0,1.,7.])
CONTROL_POINTS[7,:]=np.array([ 4.2,10,1,1.,16.])
CONTROL_POINTS[4,:]=np.array([ 5.5,5.5,0,1.,4.])
CONTROL_POINTS[13,:]=np.array([ 5.5,5.5,0.1,1.,13.])
print(CONTROL_POINTS)
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


fixed_dof=np.array([0,9])
fixed_index=np.sort(np.concatenate((fixed_dof*3,3*fixed_dof+1,fixed_dof*3+2)))
fixed_dofs1=np.array([3,6,12,15])
fixed_index=np.sort(np.concatenate((fixed_dofs1*3,fixed_index)))

load_dof=np.array([2,8,11,17])
load_index=np.sort(load_dof*3+1)
F_E[load_index]=2.5
load_dof=np.array([5,14])
load_index=np.sort(load_dof*3+1)
F_E[load_index]=5
reduced_k=np.delete(np.delete(K_G, fixed_index, 0),fixed_index , 1)
reduced_f=np.delete(F_E, fixed_index, 0)
print(reduced_f)
U_red = np.linalg.solve(reduced_k, reduced_f)
print(U_red)
''
U_exact=[]
for i  in range(len(CONTROL_POINTS[:,0])):
    x1=CONTROL_POINTS[i,0]
    y1=CONTROL_POINTS[i,1]
    z1=CONTROL_POINTS[i,2]
    Ux= 0.009375*x1 
    Uy = -0.003125*y1
    Uz=0
    U_exact.append(Ux)
    U_exact.append(Uy)
    U_exact.append(Uz)
U_exact=np.array(U_exact).reshape(len(CONTROL_POINTS), 3)    
print(U_exact)
'''

reduced_F=np.zeros(3)
U_red=np.zeros(3)
U_red = np.linalg.solve(reduced_k, reduced_F)
n=0
#for j in loading_index:
#    U = np.insert(U, j, U_red[n])
#    n+=1
print(U)
load_index=np.array([4])
loading_index=np.sort(np.concatenate((load_index*3,3*load_index+1,load_index*3+2)))
F_E[loading_index]=-1
load_index=np.array([14])
loading_index=np.sort(np.concatenate((load_index*3,3*load_index+1,load_index*3+2)))
F_E[loading_index]=-1
red=F_E-np.matmul(K_G,U)

strain_energy_each_el=np.zeros(nel)
for i in range(0, nel):
    el_in = element_indicies[i, :]
    index=np.sort(np.concatenate((el_in*3,3*el_in+1,el_in*3+2)))
    F_element=F[index]
    U_element=U[index]
    strain_energy_each_el[i]=0.5*np.dot(F_element,U_element)
print(strain_energy_each_el)
















length=1
height=1
width=1
nx=2
ny=2
nz=2
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

CONTROL_POINTS=np.array([[0,0,0,1,0],
                        [1,0,0,1,1],
                        [0,1,0,1,2],
                        [1,1,0,1,3],
                        [0.3,0.4,0.3,1,4],#4
                        [0.6,0.3,0.3,1,5],#1
                        [0.3,0.7,0.4,1,6],#3
                        [0.7,0.7,0.4,1,7],#2
                        [0.3,0.4,0.7,1,8],#8
                        [0.7,0.3,0.6,1,9],#5
                        [0.4,0.7,0.7,1,10],#7
                        [0.6,0.6,0.6,1,11],#6
                        [0,0,1,1,12],
                        [1,0,1,1,13],
                        [0,1,1,1,14],
                        [1,1,1,1,15]])


element_index=np.array([[0,1,2,3,4,5,6,7],
                        [4,5,6,7,8,9,10,11],
                        [8,9,10,11,12,13,14,15],
                        [0,4,2,6,12,8,14,10],
                        [5,1,7,3,9,13,11,15]])
print(CONTROL_POINTS[:,:-2])

K_G=np.zeros(16*3)
F=np.zeros(16*3)
nel=5
XI_KNOTVECTOR = C.xi_knotvector()
ETA_KNOTVECTOR = C.eta_knotvector()
NETA_KNOTVECTOR = C.neta_knotvector()
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

U=[]
for i  in range(len(CONTROL_POINTS[:,0])):
    x1=CONTROL_POINTS[i,0]
    y1=CONTROL_POINTS[i,1]
    z1=CONTROL_POINTS[i,2]
    Ux=1e-3*(2*x1+y1+z1)/2
    Uy=1e-3*(x1+2*y1+z1)/2
    Uz=1e-3*(x1+y1+2*z1)/2
    #Ux=((1/100)*(3*x1+2*y1+z1))
    #Uy=((1/100)*(x1+2*y1-z1))
    #Uz=((1/100)*(-2*x1+y1+4*z1))
    U.append(Ux)
    U.append(Uy)
    U.append(Uz)
#U=np.concatenate(Ux.T,np.concatenate(Uy,Uz))
U[(13*3)]=U[(13*3+1)]=U[(13*3+2)]=0
U=np.array(U)
#CONTROL_POINTS[13,:]=np.array([ 0.6,0.4,0.3,1.,13.])
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
print('span_index',span_index)
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

dof=3
#loading_nodes=np.array([index for index,j in enumerate(CONTROL_POINTS) if CONTROL_POINTS[index,0]==0 and CONTROL_POINTS[index,0]==length])
#loading_index=np.sort(np.concatenate((loading_nodes*dof,dof*loading_nodes+1,loading_nodes*dof+2)))
#F_E[loading_index]=1
#For 3d unconstrained structure
#3 translation and 3 rotation should be present
# total of 6 zero eigen values should be present for stiffness 
U = np.linalg.solve(reduced_k, F)
#print(F.reshape(len(CONTROL_POINTS), 3))
#print(element_indicies)
strain=np.zeros(nel)
for i in range(0, nel):
    el_in = element_indicies[i, :]
    index=np.sort(np.concatenate((el_in*3,3*el_in+1,el_in*3+2)))
    F_element=F[index]
    U_element=U[index]
    strain[i]=0.5*np.dot(F_element,U_element)

print(strain)
'''