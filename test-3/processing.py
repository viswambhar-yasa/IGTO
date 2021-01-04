from Inputs import *
from geometry import knot_connectivity,controlpointassembly
from element_routine import assemble,element_routine,apply_BC 
from boundary_conditions import BC_Switcher
#import pyvtk



import numpy as np

length=48
height=12
width=1
option=3
nx=16
ny=9
nz=3
density=7850
volume_frac=0.2
pmax=5
gmax=1
rmin=1.5
load=-20

XI_DEGREE=1
ETA_DEGREE=1
NETA_DEGREE=1
    
Youngs_modulus=100000
poission_ratio=0.3


N=nx
P=ny
Q=nz


C=Inputs(length,height,width,N,P,Q,XI_DEGREE,ETA_DEGREE,NETA_DEGREE)

CONTROL_POINTS=C.crtpts_coordinates()
#print(CONTROL_POINTS)
WEIGHTS=CONTROL_POINTS[:,-1]


XI_KNOTVECTOR=C.xi_knotvector()
ETA_KNOTVECTOR=C.eta_knotvector()
NETA_KNOTVECTOR=C.neta_knotvector()
XI_SPAN,XI_KNOTCONNECTIVITY,XI_UNIKNOTS,nU=C.xi_knotspan()
ETA_SPAN,ETA_KNOTCONNECTIVITY,ETA_UNIKNOTS,nV=C.eta_knotspan()
NETA_SPAN,NETA_KNOTCONNECTIVITY,NETA_UNIKNOTS,nW=C.neta_knotspan()

print('\n')

ncp=N*P*Q
dof=3
dofcp=ncp*dof
nel=nU*nV*nW
width1=120
print('#'*width1)
fmt='{:^'+str(width1)+'}'
print(fmt.format('Dimensions of the structure \n'))
print('                   Length       :',length,'               Height     :',height,'               Width        :',width,'\n')

print('                   Xi degree    :',XI_DEGREE,'                Eta degree :',ETA_DEGREE,'                Neta degree  :',NETA_DEGREE,'\n')
print('                   NX           :',N-XI_DEGREE,'                NY         :',P-ETA_DEGREE,'                NZ           :',Q-NETA_DEGREE,'\n')
print('Number of degrees of freedom :',dofcp,'\n')

print('Number of Elements:',nel,'\n')

print('No of control points in each element:',(XI_DEGREE+1)*(ETA_DEGREE+1)*(NETA_DEGREE+1),'\n')

print('>'*width1)
print('Length of the knot vector in respective direction \n')
print('XI Vector        :',list(XI_KNOTVECTOR),'\nETA vector       :',list(ETA_KNOTVECTOR),'\nNETA vector      :', list(NETA_KNOTVECTOR),'\n')
print('<'*width1)
K_G=np.zeros((dofcp,dofcp))
F_E=np.zeros(dofcp)
U=np.zeros((dofcp))
print('#'*width1)
fmt='{:^'+str(width1)+'}'
print(fmt.format('Program has started \n'))

K_disp=True
element_indicies=controlpointassembly(N,P,Q,nU,nV,nW,XI_DEGREE,ETA_DEGREE,NETA_DEGREE,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)
span_index=knot_connectivity(N,P,Q,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)
#print(element_indicies)
print('$'*width1)
fmt='{:^'+str(width1)+'}'
print(fmt.format('Finite Element Analysis based on ISO-Geometric analysis(NURBS)\n'))
for i in range(0,nel):
    el_in=element_indicies[i,:]
    sp_in=span_index[i,:]
    X=CONTROL_POINTS[el_in,0]
    Y=CONTROL_POINTS[el_in,1]
    Z=CONTROL_POINTS[el_in,2]
    weights=CONTROL_POINTS[el_in,3]
    Uspan=XI_SPAN[sp_in[0],:]
    Vspan=ETA_SPAN[sp_in[1],:]
    Wspan=NETA_SPAN[sp_in[2],:]

    K_E,NURBS,R=element_routine(X,Y,Z,weights,Youngs_modulus,poission_ratio,Uspan,Vspan,Wspan,XI_DEGREE,XI_KNOTVECTOR,ETA_DEGREE,ETA_KNOTVECTOR,NETA_DEGREE,NETA_KNOTVECTOR)
    K_G=assemble(K_G,K_E,el_in,ncp,K_disp)
    K_disp=False


bc_disp=False
abc_disp=True
BC=BC_Switcher(CONTROL_POINTS,length,height,width,bc_disp)
fixed_dof,load_dof,fixed_pts,load_pts=BC.indirect(option)
#print(load_pts)
reduced_k,reduced_F=apply_BC(K_G,F_E,fixed_dof,load_dof,load,abc_disp)
U=np.linalg.solve(reduced_k,reduced_F)
print('Calculating Displacements \n')
#print(load_pts)
#print(np.round((reduced_k@U),10)-reduced_F)
for j in fixed_dof: 
    U=np.insert(U,j,0)
   
F_E[load_dof]=load
#print(U)
#print(F_E)
U_new=np.array((U.reshape(len(CONTROL_POINTS),3)),dtype='float64')
print('Displacements due to loading :\n')
#print(U_new)
New_control_points=CONTROL_POINTS[:,:-2]+U.reshape(len(CONTROL_POINTS),3)
Ux=U_new[:,0]
Uy=U_new[:,1]
Uz=U_new[:,2]
#print(Uy)
if option==3:
    from analytical_solution import exact_displacements
    ex_dis=np.array(exact_displacements(load,length,width,height,Youngs_modulus,poission_ratio,nx,ny))
    print('Maximum defelection of Cantilever beam-(Analytical solution)        : ',np.max(ex_dis[:,1:]))
    print('Maximum defelection of Cantilever beam-(Iso-geometric FEM solution) : ',np.max(abs(U_new)))
energy_stored=np.dot(0.5,U@F_E)
print('\nThe structure is not optimised and has 100% volume \n')
print('The strain energy of the structure               :',energy_stored)
print('$'*width1)


#print(CONTROL_POINTS,'\n')
#print(New_control_points)
#print(np.max(Uy))
'''
I=((width*(height**3))/12)
maximum_deflection=(10*(length**3))/(3*Youngs_modulus*I)
print(maximum_deflection)
'''
'''
exact_disp=exact_displacements(10,length,width,height,Youngs_modulus,nx)
#print(exact_disp)

nodes=np.array([index for index,j in enumerate(CONTROL_POINTS) if CONTROL_POINTS[index,1]==0 or CONTROL_POINTS[index,2]==0])
node_ind=np.sort(np.concatenate((dof*nodes+1,nodes*dof+2)))
#print(U[node_ind])




#x1,x2,x3=plot3d(CONTROL_POINTS)
#y1,y2,y3=plot3d(New_control_points)

fig = plt.figure()
ax = plt.axes(projection='3d')
#ax.scatter3D(x1,x2,x3, 'gray')
#ax.scatter3D(y1,y2,y3, 'yellow')
ax.scatter3D(CONTROL_POINTS[:,0],CONTROL_POINTS[:,1], CONTROL_POINTS[:,2])
ax.scatter3D(New_control_points[:,0],New_control_points[:,1], New_control_points[:,2])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
#plt.show()

x=np.array(CONTROL_POINTS[:,0])
y=np.array(CONTROL_POINTS[:,1])
z=np.array(CONTROL_POINTS[:,2])
ndim = 3
#pointsToVTKAsTIN("./open_pit", x, y, z)
#gridToVTK("./structured",x,y,z)

from mayavi import mlab
s = mlab.mesh(CONTROL_POINTS[:,0], CONTROL_POINTS[:,1], CONTROL_POINTS[:,2])
mlab.show()
'''
