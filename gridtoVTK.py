
from test_preprocessing import Inputs
import numpy as np

def VTK(CONTROL_POINTS,element_density,nU,nV,nW,file_name):
    width1=120
    print('%' * width1)
    print('\n Creating VTK file ')

    from pyevtk.hl import gridToVTK
    from numpy import array,unique,reshape

    x = array(unique(CONTROL_POINTS[:, 0]), dtype='float64')
    y = array(unique(CONTROL_POINTS[:, 1]), dtype='float64')
    z = array(unique(CONTROL_POINTS[:, 2]), dtype='float64')

    volume = (element_density.reshape((nU, nV, nW), order='F'))

    print(volume)
    path='./'+file_name
    #print('"'+path+'"')
    gridToVTK("./cantilever_beam", x, y, z, cellData={"Volume1": volume})

    print('\n VTK file generated \n')
    print('%' * width1)
    pass

