#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#----------------------------------------------------------------------#
#GRID TO VTK - Python file used to generate VTK file for visualization.
#----------------------------------------------------------------------#

from Preprocessing import Inputs
import numpy as np

def VTK(CONTROL_POINTS,element_density,nx,ny,nz,file_name):
    '''
    For visualization of the mesh, we generate a VTK file containg all the information 

    Parameters
    ----------
    CONTROL_POINTS : array
                        An array of co-ordinates of control points.

    element_density : array
                        A 1D array conatining the element density obtained after topology optimization.

    nx,ny,nz  : int
                 Number of element along x,y,z direction

    file_name : str
                    path and name of the output VTK file

    Returns
    -------
            A VTK file is generated with given file name at the given path


    '''
    width1=120
    print('%' * width1)
    print('\n Creating VTK file ')

    from pyevtk.hl import gridToVTK
    from numpy import array,unique,reshape

    x = array(unique(CONTROL_POINTS[:, 0]), dtype='float64')
    y = array(unique(CONTROL_POINTS[:, 1]), dtype='float64')
    z = array(unique(CONTROL_POINTS[:, 2]), dtype='float64')

    volume = (element_density.reshape((nx, ny, nz), order='F'))

    print(volume)
    path='./'+file_name
    #print('"'+path+'"')
    gridToVTK("./cantilever_beam", x, y, z, cellData={"Volume1": volume})

    print('\n VTK file generated \n')
    print('%' * width1)
    pass

