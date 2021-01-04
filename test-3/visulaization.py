#from geometry import knot_index,bspline_basis,controlpointassembly,knot_connectivity
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


ncp=N*P*Q
dof=3
dofcp=ncp*dof
nel=nU*nV*nW
#CONTROL_POINTS=cp(nx-1,ny-1,nz-1)
#element_indicies=controlpointassembly(N,P,Q,nU,nV,nW,XI_DEGREE,ETA_DEGREE,NETA_DEGREE,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)
#span_index=knot_connectivity(N,P,Q,XI_KNOTCONNECTIVITY,ETA_KNOTCONNECTIVITY,NETA_KNOTCONNECTIVITY)
#print(CONTROL_POINTS)
#print(element_indicies)
'''
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

import pyvista as pv
from pyvista import examples

mesh = examples.load_hexbeam()
bcpos = [(6.20, 3.00, 7.50),
         (0.16, 0.13, 2.65),
         (-0.28, 0.94, -0.21)]

p = pv.Plotter()
p.add_mesh(mesh, show_edges=True, color='white')
p.add_mesh(pv.PolyData(mesh.points), color='red',
       point_size=10, render_points_as_spheres=True)
p.camera_position = bcpos
p.show(screenshot='beam_nodes.png')

import vtk

x = [
    -1.22396, -1.17188, -1.11979, -1.06771, -1.01562, -0.963542,
    -0.911458, -0.859375, -0.807292, -0.755208, -0.703125, -0.651042,
    -0.598958, -0.546875, -0.494792, -0.442708, -0.390625, -0.338542,
    -0.286458, -0.234375, -0.182292, -0.130209, -0.078125, -0.026042,
    0.0260415, 0.078125, 0.130208, 0.182291, 0.234375, 0.286458,
    0.338542, 0.390625, 0.442708, 0.494792, 0.546875, 0.598958,
    0.651042, 0.703125, 0.755208, 0.807292, 0.859375, 0.911458,
    0.963542, 1.01562, 1.06771, 1.11979, 1.17188]

y = [
    -1.25, -1.17188, -1.09375, -1.01562, -0.9375, -0.859375,
    -0.78125, -0.703125, -0.625, -0.546875, -0.46875, -0.390625,
    -0.3125, -0.234375, -0.15625, -0.078125, 0, 0.078125,
    0.15625, 0.234375, 0.3125, 0.390625, 0.46875, 0.546875,
    0.625, 0.703125, 0.78125, 0.859375, 0.9375, 1.01562,
    1.09375, 1.17188, 1.25]

z = [
    0, 0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.75, 0.8, 0.9, 1,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
    1.7, 1.75, 1.8, 1.9, 2, 2.1,
    2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
    2.75, 2.8, 2.9, 3, 3.1, 3.2,
    3.3, 3.4, 3.5, 3.6, 3.7, 3.75,
    3.8, 3.9]

# Create a rectilinear grid by defining three arrays specifying the
# coordinates in the x-y-z directions.
xCoords = vtk.vtkFloatArray()
for i in x:
    xCoords.InsertNextValue(i)

yCoords = vtk.vtkFloatArray()
for i in y:
    yCoords.InsertNextValue(i)
    
zCoords = vtk.vtkFloatArray()
for i in z:
    zCoords.InsertNextValue(i)

# The coordinates are assigned to the rectilinear grid. Make sure that
# the number of values in each of the XCoordinates, YCoordinates,
# and ZCoordinates is equal to what is defined in SetDimensions().
#
rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions(len(x), len(y), len(z))
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

# Extract a plane from the grid to see what we've got.
plane = vtk.vtkRectilinearGridGeometryFilter()
plane.SetInputData(rgrid)
plane.SetExtent(0, 46, 16, 16, 0, 43)

rgridMapper = vtk.vtkPolyDataMapper()
rgridMapper.SetInputConnection(plane.GetOutputPort())

wireActor = vtk.vtkActor()
wireActor.SetMapper(rgridMapper)
wireActor.GetProperty().SetRepresentationToWireframe()
wireActor.GetProperty().SetColor(0, 0, 0)

# Create the usual rendering stuff.
renderer = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(renderer)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

renderer.AddActor(wireActor)
renderer.SetBackground(1, 1, 1)
renderer.ResetCamera()
renderer.GetActiveCamera().Elevation(60.0)
renderer.GetActiveCamera().Azimuth(30.0)
renderer.GetActiveCamera().Zoom(1.0)

renWin.SetSize(300, 300)

# interact with data
renWin.Render()
iren.Start()
'''
import numpy as np

x=np.array([1,1,1,1])
print(np.exp(x)+np.exp(2))

import pyvista as pv

cpos = [(-0, -0, 0), (74.8305, 89.2905, 100.0), (0.23, 0.072, 0.97)]
# read the data
grid = pv.read('cantilever_beam.vtr')
pv.set_plot_theme("night")
'''
p = pv.Plotter()
p.add_mesh(grid,opacity="geom")
p.show_grid()
p.show(auto_close=False)
viewup = [0.5, 0.5, 1]
path = p.generate_orbital_path(factor=2.0, shift=10000, viewup=viewup, n_points=36)
p.open_gif("orbit.gif")
p.orbit_on_path(path,viewup=[0, 0, 1])
p.close()
'''
#print(grid)
#pv.set_plot_theme("ParaView")


#grid.add_volume(cells, cmap="viridis", opacity=volume)
#grid.show()
# plot the data with an automatically created Plotter
#grid.plot(screenshot='cantilever_beam.png',show_edges=True,opacity="geom",show_scalar_bar=True, show_axes=True)

plotter = pv.Plotter(shape=(1, 2))
plotter.subplot(0, 0)
grid= pv.read('cantilever_beam.vtr')
#pv.set_plot_theme("night")
plotter.add_text("With opacity", font_size=10)
plotter.add_mesh(grid,opacity="linear")
plotter.show_bounds(all_edges=True)
#plotter.add_mesh(grid,opacity="linear",show_edges=True)
plotter.add_axes(interactive=True)

plotter.subplot(0, 1)
grid= pv.read('cantilever_beam.vtr')
#plotter.set_plot_theme("night")
plotter.add_text("Without opacity", font_size=10)
plotter.add_mesh(grid)
plotter.add_axes(interactive=True)
plotter.show_bounds(all_edges=True)
plotter.show(screenshot='canti.png')
plotter.show_bounds(all_edges=True)

#plotter.add_mesh(mesh, color="orange")
#plotter.show()
'''
plotter.subplot(0, 1)
plotter.add_text("Render Window 1", font_size=30)
plotter.add_mesh(pv.Cube(), show_edges=True, color="tan")
#grid.plot(screenshot='cantilever_beam.png',opacity="linear",show_scalar_bar=True, show_axes=True)
'''
'''
import numpy as np

filename = "sphere-shrinking.mp4"

mesh = pv.read('cantilever_beam.vtr')
print(mesh)
print(mesh.n_cells)
mesh.cell_arrays["data"] = np.random.random(mesh.n_cells)

plotter = pv.Plotter()
# Open a movie file
plotter.open_movie(filename)

# Add initial mesh
plotter.add_mesh(mesh, scalars="data", clim=[0, 1])
# Add outline for shrinking reference
plotter.add_mesh(mesh.outline_corners())

print('Orient the view, then press "q" to close window and produce movie')

# Render and do NOT close
plotter.show(auto_close=False)

# Run through each frame
plotter.write_frame()  # write initial data

# Update scalars on each frame
for i in range(100):
    random_points = np.random.random(mesh.points.shape)
    mesh.points = random_points * 0.01 + mesh.points * 0.99
    mesh.points -= mesh.points.mean(0)
    mesh.cell_arrays["data"] = np.random.random(mesh.n_cells)
    plotter.write_frame()  # Write this frame

# Be sure to close the plotter when finished
plotter.close()
'''