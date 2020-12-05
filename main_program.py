from geometry import trilinear_der,knot_index,bspline_basis,derbspline_basis
import numpy as np
from Inputs import *

    


def elementorder(numx,numy,numz):
    index=0
    el_order=np.zeros((numx,numz,numy))
    for i in range(numz):
        for j in range(numy):
            for k in range(numx):
                el_order[k,i,j]=index
                index+=1
    return el_order



print(XI_DEGREE,ETA_DEGREE,NETA_DEGREE)
print(XI_KNOTVECTOR)
print(ETA_KNOTVECTOR)
print(NETA_KNOTVECTOR)
WEIGHTS=np.ones(8)
DR_DX,DR_DY,DR_DZ,R =trilinear_der(0.4,0.5,0.3,WEIGHTS)
