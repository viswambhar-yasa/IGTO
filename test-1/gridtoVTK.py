'''
def gridToVTK(path, x, y, z, cellData=None, pointData=None, fieldData=None):
    """
    Write data values as a rectilinear or rectangular grid.
    Parameters
    ----------
    path : str
        name of the file without extension where data should be saved.
    x : array-like
        x coordinate axis.
    y : array-like
        y coordinate axis.
    z : array-like
        z coordinate axis.
    cellData : dict, optional
        dictionary containing arrays with cell centered data.
        Keys should be the names of the data arrays.
        Arrays must have the same dimensions in all directions and must contain
        only scalar data.
    pointData : dict, optional
        dictionary containing arrays with node centered data.
        Keys should be the names of the data arrays.
        Arrays must have same dimension in each direction and
        they should be equal to the dimensions of the cell data plus one and
        must contain only scalar data.
    fieldData : dict, optional
        dictionary with variables associated with the field.
        Keys should be the names of the variable stored in each array.
    Returns
    -------
    str
        Full path to saved file.
    Notes
    -----
    coordinates of the nodes of the grid. They can be 1D or 3D depending if
    the grid should be saved as a rectilinear or logically structured grid,
    respectively.
    Arrays should contain coordinates of the nodes of the grid.
    If arrays are 1D, then the grid should be Cartesian,
    i.e. faces in all cells are orthogonal.
    If arrays are 3D, then the grid should be logically structured
    with hexahedral cells.
    In both cases the arrays dimensions should be
    equal to the number of nodes of the grid.
    """
    # Extract dimensions
    start = (0, 0, 0)
    nx = ny = nz = 0

    if x.ndim == 1 and y.ndim == 1 and z.ndim == 1:
        nx, ny, nz = x.size - 1, y.size - 1, z.size - 1
        isRect = True
        ftype = VtkRectilinearGrid
    elif x.ndim == 3 and y.ndim == 3 and z.ndim == 3:
        s = x.shape
        nx, ny, nz = s[0] - 1, s[1] - 1, s[2] - 1
        isRect = False
        ftype = VtkStructuredGrid
    else:
        assert False
    end = (nx, ny, nz)

    w = VtkFile(path, ftype)
    w.openGrid(start=start, end=end)
    w.openPiece(start=start, end=end)

    if isRect:
        w.openElement("Coordinates")
        w.addData("x_coordinates", x)
        w.addData("y_coordinates", y)
        w.addData("z_coordinates", z)
        w.closeElement("Coordinates")
    else:
        w.openElement("Points")
        w.addData("points", (x, y, z))
        w.closeElement("Points")

    _addDataToFile(w, cellData, pointData, fieldData)
    w.closePiece()
    w.closeGrid()
    # Write coordinates
    if isRect:
        w.appendData(x).appendData(y).appendData(z)
    else:
        w.appendData((x, y, z))
    # Write data
    _appendDataToFile(w, cellData, pointData, fieldData)
    w.save()
    return w.getFileName()

'''

#!/usr/bin/env python

"""
This example shows how to create a rectilinear grid.
"""

import vtk 

def main():
    colors = vtk.vtkNamedColors()

    x = [-1.22396, -1.17188, -1.11979, -1.06771, -1.01562, -0.963542, -0.911458, -0.859375, -0.807292, -0.755208,
         -0.703125, -0.651042, -0.598958, -0.546875, -0.494792, -0.442708, -0.390625, -0.338542, -0.286458, -0.234375,
         -0.182292, -0.130209, -0.078125, -0.026042, 0.0260415, 0.078125, 0.130208, 0.182291, 0.234375, 0.286458,
         0.338542, 0.390625, 0.442708, 0.494792, 0.546875, 0.598958, 0.651042, 0.703125, 0.755208, 0.807292, 0.859375,
         0.911458, 0.963542, 1.01562, 1.06771, 1.11979, 1.17188]
    y = [-1.25, -1.17188, -1.09375, -1.01562, -0.9375, -0.859375, -0.78125, -0.703125, -0.625, -0.546875, -0.46875,
         -0.390625, -0.3125, -0.234375, -0.15625, -0.078125, 0, 0.078125, 0.15625, 0.234375, 0.3125, 0.390625, 0.46875,
         0.546875, 0.625, 0.703125, 0.78125, 0.859375, 0.9375, 1.01562, 1.09375, 1.17188, 1.25]
    z = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.75, 1.8, 1.9, 2,
         2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.75, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.75, 3.8, 3.9]
    print(len(x), len(y), len(z))

    # Create a rectilinear grid by defining three arrays specifying the
    # coordinates in the x-y-z directions.
    xCoords = vtk.vtkDoubleArray()
    for i in range(0, len(x)):
        xCoords.InsertNextValue(x[i])
    yCoords = vtk.vtkDoubleArray()
    for i in range(0, len(y)):
        yCoords.InsertNextValue(y[i])
    zCoords = vtk.vtkDoubleArray()
    for i in range(0, len(z)):
        zCoords.InsertNextValue(z[i])

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
    plane.SetExtent(0, len(x) - 1, 16, 16, 0, len(z) - 1)

    rgridMapper = vtk.vtkPolyDataMapper()
    rgridMapper.SetInputConnection(plane.GetOutputPort())

    wireActor = vtk.vtkActor()
    wireActor.SetMapper(rgridMapper)
    wireActor.GetProperty().SetColor(colors.GetColor3d("Banana"))
    wireActor.GetProperty().EdgeVisibilityOn()

    # Create the usual rendering stuff.
    renderer = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    renderer.AddActor(wireActor)
    renderer.SetBackground(colors.GetColor3d("Beige"))
    renderer.ResetCamera()
    renderer.GetActiveCamera().Elevation(60.0)
    renderer.GetActiveCamera().Azimuth(30.0)
    renderer.GetActiveCamera().Zoom(1.0)

    renWin.SetSize(640, 480)
    renWin.SetWindowName('RGrid')

    # Interact with the data.
    renWin.Render()
    iren.Start()


if __name__ == "__main__":
    main()