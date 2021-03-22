#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#-------------------------------------------------------------------------#
#VISUALIZATION - Python file built using pyvista library to plot VTK files
#-------------------------------------------------------------------------#


from Preprocessing import *
import pyvista as pv
import os

def Folder(path):
    '''
    Folder is created based on the path provided
    '''
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:
        print ('Error: Creating directory. ' +  path)


def element_density_slided1(i,CP,nx,ny,nz,element_density,optimizer):
    '''
    Function is built on pyvista to visulaize VTK files

    Parameters
    ----------
    i : int
        Loop number used in saving the plot name.
    CP : array
        Control point co-ordinates.
    nx,ny,nz : int
                No of elements along x,y,z direction.
    element_density : array
                        An 1D array containing the volume of the optimized structure.
    optimizer : str
                    Name of the optimizer used for Topology optimization.

    Returns
    -------
    A 3D plot with XY ,ZX, YZ projection are plotted and save in results folder based on loop number and optimizer.

    '''
    pv.set_plot_theme("paraview")
    
    Folder('./results/')
    path1="./results/"+optimizer+"_IGTO/"
    Folder(path1)
    path=path1+optimizer+"_cantilever_"+str(i)+".jpg"
    mesh = pv.StructuredGrid()
# Set the coordinates from the numpy array
    mesh.points = CP
# set the dimensions
    mesh.dimensions = [nx, ny, nz]
    mesh.cell_arrays["volume"] =element_density

    p = pv.Plotter(shape=(2,2))
# XYZ - show 3D scene first
    p.subplot(1,1)
    no=max(nx,ny,nz)*10
    slices=mesh.slice_along_axis(n=no, axis='x')
    p.add_mesh(slices)
    p.add_text("ISOMETRIC VIEW", font_size=10)
    p.add_mesh(slices,opacity="linear",style='surface')
    p.add_axes(interactive=True)
    p.show_grid()
    #p.camera_position = cpos
# XY
    p.subplot(0,0)
    no=ny*10
    slices=mesh.slice_along_axis(n=no, axis='y')
    p.add_mesh(slices)
    p.add_text("XY", font_size=10)
    p.add_mesh(slices,opacity="linear",style='wireframe')
    p.show_grid()
    p.add_axes(interactive=True)
    p.camera_position = 'xy'
    p.enable_parallel_projection()
# ZY
    p.subplot(0,1)
    no=ny*10
    slices=mesh.slice_along_axis(n=no, axis='y')
    p.add_mesh(slices)
    p.add_text("ZY", font_size=10)
    p.add_mesh(slices,opacity="linear",style='wireframe')
    p.add_axes(interactive=True)
    p.show_grid()
    p.camera_position = 'zy'
    p.enable_parallel_projection()
# XZ
    p.subplot(1,0)
    no=nz*10
    slices=mesh.slice_along_axis(n=no, axis='z')
    p.add_mesh(slices)
    p.add_text("XZ", font_size=10)
    p.add_mesh(slices,opacity="linear",style='wireframe')
    p.add_axes(interactive=True)
    p.show_grid()
    p.camera_position = 'xz'
    p.enable_parallel_projection()
#plotting and saving the image 
    p.show(title='3D View',screenshot=path,auto_close=True)
    pass

def deformation_plot(loop,CP,U,nx,ny,nz,element_density,optimizer):
    '''
    Function is built on pyvista to visulaize VTK files

    Parameters
    ----------
    loop : int
        Loop number used in saving the plot name.
    CP : array
        Control point co-ordinates.

    U : array
            Displacements obtained due to loading

    nx,ny,nz : int
                No of elements along x,y,z direction.
    element_density : array
                        An 1D array containing the volume of the optimized structure.
    optimizer : str
                    Name of the optimizer used for Topology optimization.

    Returns
    -------
    A XY projection is plotted and save in results folder based on loop number and optimizer.

    '''
    pv.set_plot_theme("paraview")
    plotter = pv.Plotter(shape=(1, 2))
    plotter.subplot(0, 0)
    mesh = pv.StructuredGrid()
    
    mesh.dimensions = [nx, ny, nz]
# Set the coordinates from the numpy array
    mesh.points = CP
    mesh.cell_arrays["volume"] =element_density
# set the dimensions
    no=max(nx,ny,nz)*10
    slices=mesh.slice_along_axis(n=no, axis='y')
    plotter.add_text("Without deformation", font_size=15)
    plotter.add_mesh(slices,opacity="linear",style='wireframe')
    #plotter.camera_position = 'xy'
    plotter.camera_position = 'xy'
    plotter.enable_parallel_projection()
    plotter.show_bounds(all_edges=True)
    plotter.add_axes(interactive=True)

    plotter.subplot(0, 1)
    mesh = pv.StructuredGrid()
    
    mesh.dimensions = [nx, ny, nz]
# Set the coordinates from the numpy array
    mesh.points = CP+U.reshape(len(CP),3)
    mesh.cell_arrays["volume"] =element_density
# set the dimensions
 
    #plotter.set_plot_theme("night")
    plotter.add_text("With deformation", font_size=15)
    
    no=max(nx,ny,nz)*20
    slices=mesh.slice_along_axis(n=no, axis='y')
    plotter.add_mesh(slices,opacity="linear",style='wireframe')
    #plotter.camera_position = 'xy'
    plotter.camera_position = 'xy'
    plotter.enable_parallel_projection()
    plotter.add_axes(interactive=True)
    plotter.show_bounds(all_edges=True)
    Folder('./results/')
    path1="./results/"+optimizer+"_IGTO/"
    Folder(path1)
    path=path1+optimizer+"_deformation_"+str(loop)+".jpg"
#plotting and saving the image 
    plotter.show(title='IGA analysis',screenshot=path,auto_close=True)
    pass


def mesh_vis(CP,U,nx,ny,nz,optimizer):
    '''
    Function is built on pyvista to visulaize VTK files

    Parameters
    ----------
    CP : array
        Control point co-ordinates.

    U : array
            Displacements obtained due to loading

    nx,ny,nz : int
                No of elements along x,y,z direction.
   
    optimizer : str
                    Name of the optimizer used for Topology optimization.

    Returns
    -------
    deformed and undeforemd structure is plotted and save in results folder based on loop number and optimizer.

    '''
    plotter = pv.Plotter(shape=(1, 2))
    plotter.subplot(0, 0)
    mesh = pv.StructuredGrid()
# Set the coordinates from the numpy array
    mesh.points = CP
# set the dimensions
    mesh.dimensions = [nx, ny, nz]

    plotter.add_text("Before loading", font_size=10)
    plotter.add_mesh(mesh,opacity="linear")
    plotter.show_bounds(all_edges=True)
    #plotter.add_mesh(grid,opacity="linear",show_edges=True)
    plotter.add_axes(interactive=True)

    plotter.subplot(0, 1)
    mesh = pv.StructuredGrid()
# Set the coordinates from the numpy array
    mesh.points = CP+U.reshape(len(CP),3)
# set the dimensions
    mesh.dimensions = [nx, ny, nz]

    #plotter.set_plot_theme("night")
    plotter.add_text("After deformation", font_size=10)
    plotter.add_mesh(mesh)
    plotter.add_axes(interactive=True)
    plotter.show_bounds(all_edges=True)
    #plotter.camera_position('xy')
    
    Folder('./results/')
    path="./results/"+optimizer+"_cantilever_beam.png"
#plotting and saving the image 
    plotter.show(title='IGA analysis',screenshot=path,auto_close=True)
    pass


def element_density_vis(i,CP,nx,ny,nz,element_density,optimizer,camera_pos='xy'):
    '''
    Function is built on pyvista to visulaize VTK files

    Parameters
    ----------
    i : int
        Loop number used in saving the plot name.
    CP : array
        Control point co-ordinates.
    nx,ny,nz : int
                No of elements along x,y,z direction.
    element_density : array
                        An 1D array containing the volume of the optimized structure.
    optimizer : str
                    Name of the optimizer used for Topology optimization.
    camera_pos :str,optional
                    Position of the camera . The default is 'xy'.

    Returns
    -------
    A 3D plot with XY ,ZX, YZ projection are plotted and save in results folder based on loop number and optimizer.

    '''
    pv.set_plot_theme("paraview")
    mesh = pv.StructuredGrid()
# Set the coordinates from the numpy array
    mesh.points = CP
# set the dimensions
    mesh.dimensions = [nx, ny, nz]
    mesh.cell_arrays["volume"] =element_density
    Folder('./results/')
    path1="./results/"+optimizer+"_IGTO/"
    Folder(path1)
    path=path1+optimizer+"_cantilever_"+str(i)+".jpg"
    #plotting and saving the image 
    mesh.plot(screenshot=path,opacity='linear',show_edges=True,show_axes=True,cpos=camera_pos,style='surface') 
    pass



def element_density_slided(i,CP,nx,ny,nz,element_density,optimizer,no,ax,camera_pos='xy',OFF_S=True,slices=True):
    '''
    Function is built on pyvista to visulaize VTK files

    Parameters
    ----------
    i : int
        Loop number used in saving the plot name.
    CP : array
        Control point co-ordinates.
    nx,ny,nz : int
                No of elements along x,y,z direction.
    element_density : array
                        An 1D array containing the volume of the optimized structure.
    optimizer : str
                    Name of the optimizer used for Topology optimization.
    camera_pos :str,optional
                    Position of the camera . The default is 'xy'.
    no : int
        Number of slices.
    ax : str
        slices along which axes as to be specified .
    OFF_S : bool, optional
        If True the plots are saved . The default is True.
    slices : bool, optional
         If true slices structure is plotted. The default is True.

    Returns
    -------
    A 3D plot with XY ,ZX, YZ projection are plotted and save in results folder based on loop number and optimizer.


    '''
    pv.set_plot_theme("paraview")
    mesh = pv.StructuredGrid()
# Set the coordinates from the numpy array
    mesh.points = CP
# set the dimensions
    mesh.dimensions = [nx, ny, nz]
    mesh.cell_arrays["volume"] =element_density
    Folder('./results/')
    path1="./results/"+optimizer+"_IGTO/"
    Folder(path1)
    path=path1+optimizer+"_cantilever_"+str(i)+".jpg"
    slices = mesh.slice_orthogonal()
    if slices:
        slices=mesh.slice_along_axis(n=no, axis=ax)
#plotting and saving the image 
    slices.plot(screenshot=path,opacity='linear',cpos=camera_pos,style='wireframe',off_screen=OFF_S)
    #mesh.plot(screenshot=path,opacity=1,show_edges=True,show_axes=True,cpos=camera_pos) 
    #mesh.plot(screenshot=path,opacity='linear',volume=True,show_edges=True,show_axes=True,cpos=camera_pos) 
    pass
