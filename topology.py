from processing import *
from optimization import Knearestneighbours,Optimality_criteria



Emin=1e-09
E0=1
K_G=np.zeros((dofcp,dofcp))
F_E=np.zeros(dofcp)
U=np.zeros((dofcp))
print('Optimization has started')
density_basis=np.zeros(nel)
element_density=np.ones(nel)*volume_frac
density_basis_dx=np.zeros(nel)
dcompliance=np.zeros(nel)
compliance=0
nfilter=int(nel*((2*(np.ceil(rmin)-1)+1)**2))
H,DH=Knearestneighbours(rmin,nU,nV,nW)
loop=0
change=1
print('Mass of the beam :',(length*width*height)*density)
g=0
while change >0.01 :
    if loop<20:
        penal=1
    else:
        penal=min(pmax,1.2*penal)
    K_G=np.zeros((dofcp,dofcp))
    F_E=np.zeros(dofcp)
    U=np.zeros((dofcp))
    dcompliance=np.zeros(nel)
    compliance=0
    node_density=np.ones((nel,(XI_DEGREE+1)*(ETA_DEGREE+1)*(NETA_DEGREE+1)))
    for h in range(nel):
        node_density[h,:]=node_density[h,:]*element_density[h]
    for i in range(0,nel):
        el_in=element_indicies[i,:]
        sp_in=span_index[i,:]
        X=CONTROL_POINTS[el_in,0]
        Y=CONTROL_POINTS[el_in,1]
        Z=CONTROL_POINTS[el_in,2]
        weights=CONTROL_POINTS[el_in,3]
        Uspan=XI_SPAN[sp_in[0],:]
        Vspan=ETA_SPAN[sp_in[1],:]
        Wspan=NETA_SPAN[sp_in[2],:]

        K_E,NURBS,R=element_routine(X,Y,Z,weights,Youngs_modulus,poission_ratio,Uspan,Vspan,Wspan,XI_DEGREE,XI_KNOTVECTOR,ETA_DEGREE,ETA_KNOTVECTOR,NETA_DEGREE,NETA_KNOTVECTOR)
        element_density[i]=np.dot(node_density[i,:],NURBS)
        density_basis[i]=Emin+(element_density[i]**penal)*(E0-Emin)
        K_E=density_basis[i]*K_E
        K_G=assemble(K_G,K_E,el_in,ncp)

    
    BC=BC_Switcher(CONTROL_POINTS,length,height,width)
    fixed_dof,load_dof,fixed_pts,load_pts=BC.indirect(option)
    reduced_k,reduced_F=apply_BC(K_G,F_E,fixed_dof,load_dof,load)

    U=np.matmul(np.linalg.inv(reduced_k),reduced_F)

    for j in fixed_dof:
        U=np.insert(U,j,0)
    
    F_E[load_dof]=load
    print(U@F_E)
    #energy_stored=np.dot(0.5,np.matmul(U.transpose(),np.matmul(U,KK_G)))
    for k in range(0,nel):
        el_in=element_indicies[k,:]
        sp_in=span_index[k,:]
        X=CONTROL_POINTS[el_in,0]
        Y=CONTROL_POINTS[el_in,1]
        Z=CONTROL_POINTS[el_in,2]
        weights=CONTROL_POINTS[el_in,3]
        Uspan=XI_SPAN[sp_in[0],:]
        Vspan=ETA_SPAN[sp_in[1],:]
        Wspan=NETA_SPAN[sp_in[2],:]
        
        K_E,NURBS,R=element_routine(X,Y,Z,weights,Youngs_modulus,poission_ratio,Uspan,Vspan,Wspan,XI_DEGREE,XI_KNOTVECTOR,ETA_DEGREE,ETA_KNOTVECTOR,NETA_DEGREE,NETA_KNOTVECTOR)
        element_density[i]=np.dot(node_density[i,:],NURBS)
        density_basis[k]=Emin+(element_density[k]**penal)*(E0-Emin)
        density_basis_dx[k]=-penal*(element_density[k]**(penal-1))*(E0-Emin)
        dof_index=np.sort(np.concatenate((el_in*dof,dof*el_in+1,dof*el_in+2)))
        #compliance+=np.transpose(U[dof_index])@(density_basis[k]*K_E)@U[dof_index]
        dcompliance[k]=np.transpose(U[dof_index])@(density_basis_dx[k]*K_E)@U[dof_index]
    compliance+=np.transpose(U)@(K_G)@U
    dv=np.ones(nel)
    dcompliance=np.matmul(H,dcompliance/DH)
    dv=np.matmul(H,dv/DH)
    old_el_density=element_density
    element_density_updated= Optimality_criteria(dcompliance,dv,H,DH,element_density,nel,volume_frac)
    element_density=element_density_updated
    change=np.linalg.norm(element_density-old_el_density)

    if loop==0:
        initial_compliance=compliance
    loop+=1
    print('  ITR: ', loop,'      OBJ: ',compliance,'      Vol: ',  ((volume_frac)),'     CHG: ',change)
print('\n Optimized mass of the beam :',(length*width*height)*density*np.mean(element_density))
energy_stored=np.dot(0.5,np.matmul(U.transpose(),np.matmul(U,K_G)))
print('\n Optimised structure Compliance',initial_compliance,'\n')
final_compliance=np.dot(U,F_E)
print('\n Optimised structure Compliance',final_compliance,'\n')
print(np.round(element_density,2))