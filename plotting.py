from Inputs import *
import pyvista as pv
import vtk
from geometry import controlpointassembly


length=48
height=12
width=6
option=3
nx=3
ny=2
nz=2
density=7850
volume_frac=0.5
pmax=3
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

WEIGHTS=CONTROL_POINTS[:,-1]


XI_KNOTVECTOR=C.xi_knotvector()
ETA_KNOTVECTOR=C.eta_knotvector()
NETA_KNOTVECTOR=C.neta_knotvector()
XI_SPAN,XI_KNOTCONNECTIVITY,XI_UNIKNOTS,nU=C.xi_knotspan()
ETA_SPAN,ETA_KNOTCONNECTIVITY,ETA_UNIKNOTS,nV=C.eta_knotspan()
NETA_SPAN,NETA_KNOTCONNECTIVITY,NETA_UNIKNOTS,nW=C.neta_knotspan()



ncp=N*P*Q
dof=3
dofcp=ncp*dof
nel=nU*nV*nW

points=CONTROL_POINTS[:,:-2]
print(points)
cell=np.zeros((nel,(XI_DEGREE+1)*(ETA_DEGREE+1)*(NETA_DEGREE+1)))
cell=controlpointassembly(N,P,Q,nU,nV,nW,XI_DEGREE,ETA_DEGREE,NETA_DEGREE,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)


print(cell)
print(cell.ravel())
celltype=np.empty(nel,dtype=np.uint8)
celltype[:] = vtk.VTK_HEXAHEDRON

grid = pv.UnstructuredGrid(cell, celltype, points)
_ = grid.plot(show_edges=True)
'''
import numpy as np
from pyvista import examples


dataset = examples.
dataset.set_active_scalars("Spatial Cell Data")

print(dataset)
'''