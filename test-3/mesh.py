from Inputs import *
import pyvista as pv
length=48
height=12
width=1
option=3
nx=3
ny=3
nz=3
density=7850
volume_frac=0.5
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
X=CONTROL_POINTS[:,0]
Y=CONTROL_POINTS[:,1]
Z=CONTROL_POINTS[:,2]
mesh = pv.StructuredGrid()
# Set the coordinates from the numpy array
mesh.points = CONTROL_POINTS[:,:-2]
# set the dimensions
mesh.dimensions = [nx, ny, nz]
mesh.cell_arrays["volume"] =[1,0,1,0,1,1,0,1]
# and then inspect it!
mesh.plot(show_edges=True, show_grid=True,opacity='linear', cpos="xy")
print(mesh)