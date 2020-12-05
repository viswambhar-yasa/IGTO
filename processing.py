from Inputs import *
from geometry import knot_connectivity,controlpointassembly
from element_routine import assemble,element_routine,apply_BC 
from boundary_conditions import BC_Switcher
#import pyvtk
from analytical_solution import exact_displacements
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


import numpy as np

length=1
height=1
width=1
option=0
'''
provide the number of elements in each direction

'''
nx=5
ny=4
nz=4

'''
Provide the  information for knot vector in each direction
'''

N=nx
P=ny
Q=nz

XI_DEGREE=1
ETA_DEGREE=1
NETA_DEGREE=1

C=Inputs(length,height,width,N,P,Q,XI_DEGREE,ETA_DEGREE,NETA_DEGREE)

CONTROL_POINTS=C.crtpts_coordinates()

WEIGHTS=CONTROL_POINTS[:,-1]


XI_KNOTVECTOR=C.xi_knotvector()
ETA_KNOTVECTOR=C.eta_knotvector()
NETA_KNOTVECTOR=C.neta_knotvector()
XI_SPAN,XI_KNOTCONNECTIVITY,XI_UNIKNOTS,nU=C.xi_knotspan()
ETA_SPAN,ETA_KNOTCONNECTIVITY,ETA_UNIKNOTS,nV=C.eta_knotspan()
NETA_SPAN,NETA_KNOTCONNECTIVITY,NETA_UNIKNOTS,nW=C.neta_knotspan()


Youngs_modulus=100000
poission_ratio=0.3

ncp=N*P*Q
dof=3
dofcp=ncp*dof
nel=nU*nV*nW
print('\n Dimensions of the structure \n')
print('Length:',length,'    Height :',height,'     Width: ',width,'\n')

print('Length of the knot vector in respective direction \n')
print('XI Vector    :',XI_KNOTVECTOR,'\nETA vector      :',ETA_KNOTVECTOR,'\nNETA vector    :', NETA_KNOTVECTOR,'\n')
print('Xi degree: ',XI_DEGREE,'     Eta degree: ',ETA_DEGREE,'      Neta degree: ',NETA_DEGREE,'\n')
print('Number of degrees of freedom :',dofcp,'\n')

print('Number of Elements:',nel,'\n')
print('Number of elements in each direction \n')
print('NX :',N-XI_DEGREE,'   NY:',P-ETA_DEGREE,'   NZ:',Q-NETA_DEGREE,'\n')
print('No of control points in each element:',(XI_DEGREE+1)*(ETA_DEGREE+1)*(NETA_DEGREE+1),'\n')


K_G=np.zeros((dofcp,dofcp))
F_E=np.zeros(dofcp)
U=np.zeros((dofcp))
print('Program has started')


element_indicies=controlpointassembly(N,P,Q,nU,nV,nW,XI_DEGREE,ETA_DEGREE,NETA_DEGREE,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)
span_index=knot_connectivity(N,P,Q,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)

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
    K_G=assemble(K_G,K_E,el_in,ncp)



BC=BC_Switcher(CONTROL_POINTS,length,height,width)
fixed_dof,load_dof,fixed_pts,load_pts=BC.indirect(option)
reduced_k,reduced_F=apply_BC(K_G,F_E,fixed_dof,load_dof,10)
U=np.matmul(np.linalg.inv(reduced_k),reduced_F)

#print(np.round((reduced_k@U),10)-reduced_F)
for j in fixed_dof: 
    U=np.insert(U,j,0)
F_E[load_dof]=P
print(U)
U_new=np.array((U.reshape(len(CONTROL_POINTS),3)),dtype='float64')
#print(U_new)
New_control_points=CONTROL_POINTS[:,:-2]+U.reshape(len(CONTROL_POINTS),3)
Uy=U_new[:,1]
print(np.max(Uy))
I=(width*(height**3))/12
maximum_deflection=(10*(length**3))/(3*Youngs_modulus*I)
print(maximum_deflection)
energy_stored=np.dot(0.5,np.matmul(U.transpose(),np.matmul(U,K_G)))
print('\n strain energy stored in the Structure',energy_stored,'\n')


exact_disp=exact_displacements(10,length,width,height,Youngs_modulus,nx)
#print(exact_disp)

nodes=np.array([index for index,j in enumerate(CONTROL_POINTS) if CONTROL_POINTS[index,1]==0 or CONTROL_POINTS[index,2]==0])
node_ind=np.sort(np.concatenate((dof*nodes+1,nodes*dof+2)))
#print(U[node_ind])
'''
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