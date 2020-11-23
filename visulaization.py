from geometry import knot_index,bspline_basis


def NURBS_volume(Ux,Uy,Uz,xdegree,xknotvector,ydegree,yknotvector,zdegree,zknotvector,weights,control_points):  
    '''
    Return a point based 

    Parameters
    ----------
    Ux : TYPE
        DESCRIPTION.
    Uy : TYPE
        DESCRIPTION.
    Uz : TYPE
        DESCRIPTION.
    xdegree : TYPE
        DESCRIPTION.
    xknotvector : TYPE
        DESCRIPTION.
    ydegree : TYPE
        DESCRIPTION.
    yknotvector : TYPE
        DESCRIPTION.
    zdegree : TYPE
        DESCRIPTION.
    zknotvector : TYPE
        DESCRIPTION.
    weights : TYPE
        DESCRIPTION.
    control_points : TYPE
        DESCRIPTION.

    Returns
    -------
    volume_points : TYPE
        DESCRIPTION.

    '''
    
    xspan=len(xknotvector)-xdegree-1
    yspan=len(yknotvector)-ydegree-1
    x_index=knot_index(xdegree,Ux,xknotvector)
    Nx=bspline_basis(x_index,xdegree,Ux,xknotvector)
    y_index=knot_index(ydegree,Uy,yknotvector)
    Ny=bspline_basis(y_index,ydegree,Uy,yknotvector)
    z_index=knot_index(zdegree,Uz,zknotvector)
    Nz=bspline_basis(z_index,zdegree,Uz,zknotvector)
    Vi=0
    Wi=0
    for k in range(zdegree+1):
        zind=z_index-ydegree+k
        for j in range(ydegree+1):
            yind=y_index-ydegree+j
            for i in range(xdegree+1):
                cind=x_index-xdegree+i+(xspan)*((yspan*zind)+yind)
                Wi=Wi+Nx[i]*Ny[j]*Nz[k]*weights[cind]
                Vi=Vi+Nx[i]*Ny[j]*Nz[k]*control_points[cind]*weights[cind]             
    volume_points=Vi/Wi
    return volume_points
