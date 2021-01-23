# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 22:24:09 2021

@author: VISWAMBHAR YASA
"""
import pyvista as pv
import vtk
import numpy as np
from Preprocessing import *
from geometry import knot_connectivity, controlpointassembly
from boundary_conditions import BC_Switcher
from visulaization import mesh_vis
import timeit

start = timeit.default_timer()

# All the program statements
stop = timeit.default_timer()
execution_time = stop - start

print("Program Executed in ",execution_time)

length=1
height=1
width=1
option=3
nx=3
ny=3
nz=3


XI_DEGREE=2
ETA_DEGREE=2
NETA_DEGREE=2
    

N = nx
P = ny
Q = nz

C = Inputs(length, height, width, N, P, Q, XI_DEGREE, ETA_DEGREE, NETA_DEGREE)

CONTROL_POINTS = C.crtpts_coordinates()
print(CONTROL_POINTS)
# print(CONTROL_POINTS)
WEIGHTS = CONTROL_POINTS[:, -1]

XI_KNOTVECTOR = C.xi_knotvector()
ETA_KNOTVECTOR = C.eta_knotvector()
NETA_KNOTVECTOR = C.neta_knotvector()
XI_SPAN, XI_KNOTCONNECTIVITY, XI_UNIKNOTS, nU = C.xi_knotspan()
ETA_SPAN, ETA_KNOTCONNECTIVITY, ETA_UNIKNOTS, nV = C.eta_knotspan()
NETA_SPAN, NETA_KNOTCONNECTIVITY, NETA_UNIKNOTS, nW = C.neta_knotspan()

element_indicies = controlpointassembly(N, P, Q, nU, nV, nW, XI_DEGREE, ETA_DEGREE, NETA_DEGREE, XI_KNOTCONNECTIVITY,
                                        ETA_KNOTCONNECTIVITY, NETA_KNOTCONNECTIVITY)

print(element_indicies)


import vtk

def main():
    colors = vtk.vtkNamedColors()

    # create polyhedron (cube)
    # The point Ids are: [0, 1, 2, 3, 4, 5, 6, 7]

    points = vtk.vtkPoints()
    points.InsertNextPoint(-1.0, -1.0, -1.0)
    points.InsertNextPoint(1.0, -1.0, -1.0)
    points.InsertNextPoint(1.0, 1.0, -1.0)
    points.InsertNextPoint(-1.0, 1.0, -1.0)
    points.InsertNextPoint(-1.0, -1.0, 1.0)
    points.InsertNextPoint(1.0, -1.0, 1.0)
    points.InsertNextPoint(1.0, 1.0, 1.0)
    points.InsertNextPoint(-1.0, 1.0, 1.0)

    # These are the point ids corresponding to each face.
    faces = [[0, 3, 2, 1], [0, 4, 7, 3], [4, 5, 6, 7], [5, 1, 2, 6], [0, 1, 5, 4], [2, 3, 7, 6]]
    faceId = vtk.vtkIdList()
    faceId.InsertNextId(6)  # Six faces make up the cell.
    for face in faces:
        faceId.InsertNextId(len(face))  # The number of points in the face.
        [faceId.InsertNextId(i) for i in face]

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(vtk.VTK_POLYHEDRON, faceId)

    # Here we write out the cube.
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName('polyhedron.vtu')
    writer.SetDataModeToAscii()
    writer.Update()

    # Create a mapper and actor
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(ugrid)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(
        colors.GetColor3d('Silver'))

    # Visualize
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName('Polyhedron')
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('Salmon'))
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(30)
    renderer.GetActiveCamera().Elevation(30)
    renderWindow.Render()
    renderWindowInteractor.Start()


if __name__ == '__main__':
    main()
    
# these points will all be shared between the cells
points = np.array([[0. , 0. , 0. ],
                   [1. , 0. , 0. ],
                   [0.5, 0. , 0. ],
                   [1. , 1. , 0. ],
                   [1. , 0.5, 0. ],
                   [0. , 1. , 0. ],
                   [0.5, 1. , 0. ],
                   [0. , 0.5, 0. ],
                   [0.5, 0.5, 0. ],
                   [1. , 0. , 0.5],
                   [1. , 0. , 1. ],
                   [0. , 0. , 0.5],
                   [0. , 0. , 1. ],
                   [0.5, 0. , 0.5],
                   [0.5, 0. , 1. ],
                   [1. , 1. , 0.5],
                   [1. , 1. , 1. ],
                   [1. , 0.5, 0.5],
                   [1. , 0.5, 1. ],
                   [0. , 1. , 0.5],
                   [0. , 1. , 1. ],
                   [0.5, 1. , 0.5],
                   [0.5, 1. , 1. ],
                   [0. , 0.5, 0.5],
                   [0. , 0.5, 1. ],
                   [0.5, 0.5, 0.5],
                   [0.5, 0.5, 1. ]])


# Each cell in the cell array needs to include the size of the cell
# and the points belonging to the cell.  In this example, there are 8
# hexahedral cells that have common points between them.
cells = np.array([[ 8,  0,  2,  8,  7, 11, 13, 25, 23],
                  [ 8,  2,  1,  4,  8, 13,  9, 17, 25],
                  [ 8,  7,  8,  6,  5, 23, 25, 21, 19],
                  [ 8,  8,  4,  3,  6, 25, 17, 15, 21],
                  [ 8, 11, 13, 25, 23, 12, 14, 26, 24],
                  [ 8, 13,  9, 17, 25, 14, 10, 18, 26],
                  [ 8, 23, 25, 21, 19, 24, 26, 22, 20],
                  [ 8, 25, 17, 15, 21, 26, 18, 16, 22]]).ravel()

# each cell is a VTK_HEXAHEDRON
celltypes = np.empty(8, dtype=np.uint8)
celltypes[:] = vtk.VTK_HEXAHEDRON
celltypes[:]   =   vtk.VTK_POLYHEDRON

# the offset array points to the start of each cell (via flat indexing)
offset = np.array([ 0, 9, 18, 27, 36, 45, 54, 63])

# Effectively, when visualizing a VTK unstructured grid, it will
# sequentially access the cell array by first looking at each index of
# cell array (based on the offset array), and then read the number of
# points based on the first value of the cell.  In this case, the
# VTK_HEXAHEDRON is described by 8 points.

# for example, the 5th cell would be accessed by vtk with:
start_of_cell = offset[4]
n_points_in_cell = cells[start_of_cell]
indices_in_cell = cells[start_of_cell + 1: start_of_cell + n_points_in_cell + 1]
print(indices_in_cell)

# if you are using VTK 9.0 or newer, you do not need to input the offset array:
# grid = pv.UnstructuredGrid(cells, celltypes, points)

# if you are not using VTK 9.0 or newer, you must use the offset array
# Alternate versions:
grid = pv.UnstructuredGrid({vtk.VTK_HEXAHEDRON: cells.reshape([-1, 9])[:, 1:]}, points)
grid = pv.UnstructuredGrid({vtk.VTK_HEXAHEDRON: np.delete(cells, np.arange(0, cells.size, 9))}, points)

# plot the grid (and suppress the camera position output)
_ = grid.plot(show_edges=True,style='wireframe')