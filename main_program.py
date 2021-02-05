import sys
import subprocess
import pkg_resources
import pip
import time 

TYELLOW =  '\033[33;1m' 
TGREEN = '\033[32;1m'
TRED='\033[31;1m'
TBLUE = '\033[34;1m'
ENDC = '\033[m' 

main_program_start=time.time()
for package in ['numpy', 'matplotlib', 'pyvista', 'pyEVTK']:
    try:
        dist = pkg_resources.get_distribution(package)
        print(TGREEN+'{} ({}) is installed'.format(dist.key, dist.version)+ENDC)
    except pkg_resources.DistributionNotFound:
        print(TRED+'{} is NOT installed'.format(package)+ENDC)
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
    from processing import *

from plotting import *
plotting(ii,CC,element_density,optimizer,option)

