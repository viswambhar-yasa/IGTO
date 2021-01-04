import numpy as np
import matplotlib.pyplot as plt
from Inputs import Inputs
from geometry import controlpointassembly
#import sys
#sys.setrecursionlimit(5000)
def knotvector(ncontrol_points,degree):
    size=ncontrol_points+degree+1
    knotvector=np.zeros(size)
    if ncontrol_points >(degree):
        for i in range(size):
            if i<=degree:
                knotvector[i]=0
            elif i>degree and i<(size-degree):
                value=i-degree
                knotvector[i]=value
            else:
                knotvector[i]=ncontrol_points-degree
    else:
        raise Exception('Require more control points to define the knotvector')
    knotvector=knotvector/knotvector[-1]
    return knotvector

def CNT(nx,ny,nz):

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

length=1
height=1
width=1

'''
provide the number of elements in each direction

'''
nx=2
ny=2
nz=2

'''
Provide the  information for knot vector in each direction
'''

N=nx
P=ny
Q=nz

XI_DEGREE=1
ETA_DEGREE=1
NETA_DEGREE=1

C=Inputs(length,height,width,nx,ny,nz,XI_DEGREE,ETA_DEGREE,NETA_DEGREE)

CONTROL_POINTS=C.crtpts_coordinates()


WEIGHTS=CONTROL_POINTS[:,-1]


XI_KNOTVECTOR=C.xi_knotvector()
ETA_KNOTVECTOR=C.eta_knotvector()
NETA_KNOTVECTOR=C.neta_knotvector()
XI_SPAN,XI_KNOTCONNECTIVITY,XI_UNIKNOTS,nU=C.xi_knotspan()
ETA_SPAN,ETA_KNOTCONNECTIVITY,ETA_UNIKNOTS,nV=C.eta_knotspan()
NETA_SPAN,NETA_KNOTCONNECTIVITY,NETA_UNIKNOTS,nW=C.neta_knotspan()

CP=CNT((nx-1),(ny-1),(nz-1))
element_indicies=controlpointassembly(N,P,Q,nU,nV,nW,XI_DEGREE,ETA_DEGREE,NETA_DEGREE,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)
print(element_indicies)
import meshio
points =  CP[:,:-2]
print(points)
#np.array(element_indicies)
cells = {
    "hexahedron": np.array([[0,1,2,3,4,5,6,7]])
    }
meshio.write_points_cells(
    "foo.vtk",
    points,
    cells
    # Optionally provide extra data on points, cells, etc.
    # point_data=point_data,
    # cell_data=cell_data,
    # field_data=field_data
)

