import numpy as np
import matplotlib.pyplot as plt
def knotvector(ncontrol_points,degree):
    size=ncontrol_points+degree+1
    knotvector=np.zeros(size)
    if ncontrol_points >(degree):
        for i in range(size):
            if i<=degree:
                knotvector[i]=0
            elif i>degree and i<(size-degree):
                value=i-degree
                knotvector[i]=value
            else:
                knotvector[i]=ncontrol_points-degree
    else:
        raise Exception('Require more control points to define the knotvector')
    knotvector=knotvector/knotvector[-1]
    return knotvector