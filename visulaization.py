from geometry import knot_index,bspline_basis,controlpointassembly,knot_connectivity
#from processing import *
#from topology import *
from Inputs import *
'''
def NURBS_volume(Ux,Uy,Uz,control_points,weights=WEIGHTS,xdegree=XI_DEGREE,xknotvector=XI_KNOTCONNECTIVITY,ydegree=ETA_DEGREE,yknotvector=ETA_KNOTVECTOR,zdegree=NETA_DEGREE,zknotvector=NETA_KNOTVECTOR):  
    
    Return a point based 

    Parameters
    ----------
    Ux : TYPE
        DESCRIPTION.
    Uy : TYPE
        DESCRIPTION.
    Uz : TYPE
        DESCRIPTION.
    xdegree : TYPE
        DESCRIPTION.
    xknotvector : TYPE
        DESCRIPTION.
    ydegree : TYPE
        DESCRIPTION.
    yknotvector : TYPE
        DESCRIPTION.
    zdegree : TYPE
        DESCRIPTION.
    zknotvector : TYPE
        DESCRIPTION.
    weights : TYPE
        DESCRIPTION.
    control_points : TYPE
        DESCRIPTION.

    Returns
    -------
    volume_points : TYPE
        DESCRIPTION.

    
    
    xspan=len(xknotvector)-xdegree-1
    yspan=len(yknotvector)-ydegree-1
    x_index=knot_index(xdegree,Ux,xknotvector)
    Nx=bspline_basis(x_index,xdegree,Ux,xknotvector)
    y_index=knot_index(ydegree,Uy,yknotvector)
    Ny=bspline_basis(y_index,ydegree,Uy,yknotvector)
    z_index=knot_index(zdegree,Uz,zknotvector)
    Nz=bspline_basis(z_index,zdegree,Uz,zknotvector)
    Vi=0
    Wi=0
    for k in range(zdegree+1):
        zind=z_index-ydegree+k
        for j in range(ydegree+1):
            yind=y_index-ydegree+j
            for i in range(xdegree+1):
                cind=x_index-xdegree+i+(xspan)*((yspan*zind)+yind)
                Wi=Wi+Nx[i]*Ny[j]*Nz[k]*weights[cind]
                Vi=Vi+Nx[i]*Ny[j]*Nz[k]*control_points[cind]*weights[cind]             
    volume_points=Vi/Wi
    return volume_points


def plot3d(CP,steps=20,xknotvector=XI_KNOTVECTOR,yknotvector=ETA_KNOTVECTOR,zknotvector=NETA_KNOTCONNECTIVITY):
    u_min = xknotvector[0]
    u_max = xknotvector[-1]
    v_min = yknotvector[0]
    v_max = yknotvector[-1]
    w_min = yknotvector[0]
    w_max = yknotvector[-1]
    #xstep_size = (u_max-u_min)/(steps-1)
    #ystep_size = (v_max-v_min)/(steps-1)
    #zstep_size = (w_max-w_min)/(steps-1)    
    a = np.zeros(steps*steps*steps)
    b = np.zeros(steps*steps*steps)
    c = np.zeros(steps*steps*steps)
    i=0
    tol=0.00001
    xn=np.linspace(u_min,u_max-tol,steps)
    yn=np.linspace(v_min,v_max-tol,steps)
    zn=np.linspace(w_min,w_max-tol,steps)
    for t in xn:
        for p in yn:
            for r in zn:
                C=NURBS_volume(t,p,r,CP)
                a[i]=C[0]
                b[i]=C[1]
                c[i]=C[2]
                i=i+1
    return a,b,c
'''


length=1
height=1
width=1
option=3
nx=5
ny=5
nz=5
density=7850
volume_frac=0.5
pmax=3
rmin=1.5
load=-20

XI_DEGREE=3
ETA_DEGREE=4
NETA_DEGREE=3
    
Youngs_modulus=100000
poission_ratio=0.3


N=nx
P=ny
Q=nz


C=Inputs(length,height,width,N,P,Q,XI_DEGREE,ETA_DEGREE,NETA_DEGREE)

CONTROL_POINTS=C.crtpts_coordinates()

WEIGHTS=CONTROL_POINTS[:,-1]


XI_KNOTVECTOR=C.xi_knotvector()
ETA_KNOTVECTOR=C.eta_knotvector()
NETA_KNOTVECTOR=C.neta_knotvector()
XI_SPAN,XI_KNOTCONNECTIVITY,XI_UNIKNOTS,nU=C.xi_knotspan()
ETA_SPAN,ETA_KNOTCONNECTIVITY,ETA_UNIKNOTS,nV=C.eta_knotspan()
NETA_SPAN,NETA_KNOTCONNECTIVITY,NETA_UNIKNOTS,nW=C.neta_knotspan()
def cp(nx,ny,nz):
    '''
     A function which return the control point co-ordinates
     '''
    beam_coordinates=[]
    index=0
    for k in range(nx+1):
        for i in range(nz+1):
            for j in range(ny+1):
                beam_coordinates.append([(k)*(length/nx),j*(height/ny),i*(width/nz),1,index])
                index+=1
    return np.array(beam_coordinates)


ncp=N*P*Q
dof=3
dofcp=ncp*dof
nel=nU*nV*nW
#CONTROL_POINTS=cp(nx-1,ny-1,nz-1)
element_indicies=controlpointassembly(N,P,Q,nU,nV,nW,XI_DEGREE,ETA_DEGREE,NETA_DEGREE,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)
span_index=knot_connectivity(N,P,Q,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)
print(CONTROL_POINTS)
print(element_indicies)
from pyevtk.hl import unstructuredGridToVTK,gridToVTK
from pyevtk.vtk import VtkQuad
import numpy as np 
import random as rnd 
#element_indicies=np.array([[0,1,3,2,7,5,4,6]])
 # Dimensions 
x=(CONTROL_POINTS[:,0]).ravel()
y=(CONTROL_POINTS[:,1]).ravel()
z=(CONTROL_POINTS[:,2]).ravel()
connectivity=element_indicies.ravel()
print(connectivity)
offset=element_indicies[:,-1:].ravel()
print(offset)
ctype = np.ones(nel)
ctype=ctype*VtkQuad.tid
#cellData=element_density.ravel()
#print(cellData)
unstructuredGridToVTK("./v13", x, y, z,connectivity,offset,ctype)

gridToVTK("./v1",x,y,z
    #cellData={"pressure": pressure},
    #pointData={"temp": temp},
)
print('VTK generated')