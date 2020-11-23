from Inputs import Inputs




if __name__ == "__main__":
    '''
    Give the required dimension of the beam 
    '''
    length=1
    height=1
    width=1

    '''
    provide the number of elements in each direction

    '''
    nx=2
    ny=2
    nz=2

    '''
    Provide the  information for knot vector in each direction
    '''
    n=nx+1
    p=ny+1
    q=nz+1

    xidegree=2
    etadegree=2
    netadegree=2

    C=Inputs(length,height,width,nx,ny,nz,xidegree,etadegree,netadegree)
    global control_points
    control_points=C.crtpts_coordinates()
    global weights
    weights=control_points[:,-1]
    global xi_knotvector
    global eta_knotvector
    global neta_knotvector

    xi_knotvector=C.xi_knotvector()
    eta_knotvector=C.eta_knotvector()
    neta_knotvector=C.neta_knotvector()
