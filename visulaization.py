from Preprocessing import *
import pyvista as pv
import os

def Folder(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:
        print ('Error: Creating directory. ' +  path)


def element_density_slided1(i,CP,nx,ny,nz,element_density,optimizer):
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

    p.show(title='3D View',screenshot=path,auto_close=True)
    pass

def deformation_plot(loop,CP,U,nx,ny,nz,element_density,optimizer):
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
    plotter.show(title='IGA analysis',screenshot=path,auto_close=True)
    pass


def mesh_vis(CP,U,nx,ny,nz,optimizer):
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
    plotter.show(title='IGA analysis',screenshot=path,auto_close=True)
    pass


def element_density_vis(i,CP,nx,ny,nz,element_density,optimizer,camera_pos='xy'):
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
    
    mesh.plot(screenshot=path,opacity='linear',show_edges=True,show_axes=True,cpos=camera_pos,style='surface') 
    pass



def element_density_slided(i,CP,nx,ny,nz,element_density,optimizer,no,ax,camera_pos='xy',OFF_S=True,slices=True):
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
    slices.plot(screenshot=path,opacity='linear',cpos=camera_pos,style='wireframe',off_screen=OFF_S)
    #mesh.plot(screenshot=path,opacity=1,show_edges=True,show_axes=True,cpos=camera_pos) 
    #mesh.plot(screenshot=path,opacity='linear',volume=True,show_edges=True,show_axes=True,cpos=camera_pos) 
    pass
