import sys
import subprocess
import pkg_resources
import pip


for package in ['numpy', 'matplotlib', 'pyvista', 'pyEVTK']:
    try:
        dist = pkg_resources.get_distribution(package)
        print('{} ({}) is installed'.format(dist.key, dist.version))
    except pkg_resources.DistributionNotFound:
        print('{} is NOT installed'.format(package))
        pip.main(['install', package])
 
import os

def Folder(path):
    '''
    

    Parameter
    ----------
    path : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:
        print('Error: Creating directory. ' + path)


Folder('./results/')

if __name__ == "__main__":
    from topology import *

from plotting import *

plotting(ii,CC,element_density,optimizer,option)

