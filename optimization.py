import numpy as np

def Knearestneighbours(rmin,nelx,nely,nelz): ##modify
    nel=nelx*nely*nelz
    H=np.zeros((nel,nel))
    for i in range(nelz):
        for j in range(nely):
            for k in range(nelx):
                r=int(i*nelx*nely+j*nelx+k)
                imin=int(max(i-(np.ceil(rmin)-1),0))
                imax=int(min(i+np.ceil(rmin),nelz))
                jmin=int(max(j-(np.ceil(rmin)-1),0))
                jmax=int(min(j+np.ceil(rmin),nely))
                kmin=int(max(k-(np.ceil(rmin)-1),0))
                kmax=int(min(k+np.ceil(rmin),nelx))
                for ii in range(imin,imax): 
                    for jj in range(jmin,jmax):
                        for kk in range(kmin,kmax):
                            c= int(ii*nelx*nely+jj*nelx+kk)
                            distance=np.sqrt((i-ii)**2+(j-jj)**2+(k-kk)**2)
                            H[r,c]=np.maximum(0,(rmin-distance))
    H=np.array(H)
    DH=np.sum(H,1)
    return H,DH  

def Optimality_criteria(dC,dV,H,DH,X,nel,vol_fra,tol=0.001,move=0.2,initial_value=0,final_value=1e09):
    convergence_criteria=(final_value-initial_value)/(final_value+initial_value)
    X_updated=np.zeros(len(X))
    X=np.array(X)
    while convergence_criteria>=tol:
        lamdba=0.5*(initial_value+final_value)
        Be=np.sqrt(-dC/(dV*lamdba))
        X_updated[:]= np.maximum(0.0,np.maximum(X-move,np.minimum(1.0,np.minimum(X+move,X*Be))))
        moving_condition=np.sum(X_updated)
        if moving_condition>vol_fra*nel:
            initial_value=lamdba
        else:
            final_value=lamdba
        convergence_criteria=(final_value-initial_value)/(final_value+initial_value)    
        
    return X_updated

