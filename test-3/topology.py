from processing import *
from optimization import Knearestneighbours,optimality_criteria,Heavyside_filter,Moving_asymptoes
import matplotlib.pyplot as plt


Emin=1e-09
E0=1
K_G=np.zeros((dofcp,dofcp))
F_E=np.zeros(dofcp)
U=np.zeros((dofcp))
print('+'*width1)
fmt='{:^'+str(width1)+'}'
print(fmt.format('Structural Topology optimization using IGA \n'))
print('\n                                            Optimization has started \n')

print('                                     Density of the material       :',density)
print('                                     Youngs Modulus                :',Youngs_modulus)
print('                                     Poission ratio                :',poission_ratio,'\n')
fmt='{:^'+str(width1)+'}'
print(fmt.format('The percentage of volume which has to remain after optimization \n'))
print('                                          ',volume_frac)
print('+'*width1)
density_basis=np.zeros(nel)
ele_den_filter=np.zeros(nel)
element_density=np.ones(nel)*volume_frac
density_basis_dx=np.zeros(nel)
dcompliance=np.zeros(nel)
compliance=0
nfilter=int(nel*((2*(np.ceil(rmin)-1)+1)**2))
H,DH=Knearestneighbours(rmin,nU,nV,nW)
loop=0
change=1
g=1
CC=[]
ii=[]
VV=[]
optimizer='MMA'
penal=max(15*((1-poission_ratio)/(7-5*poission_ratio)),(3/2)*((1-poission_ratio)/(1-2*poission_ratio)))
penal=3.5
oc_disp=True
bc_disp=True
fil_disp=True
max_iterations=250
X=[]#np.zeros(nel,max_iterations)
beta=1
filter_N=np.zeros(((XI_DEGREE+1)*(ETA_DEGREE+1)*(NETA_DEGREE+1),nel))
Xmin=np.zeros(nel)
Xmax=np.ones(nel)
Lower=Xmin
Upper=Xmin
E1=element_density
E2=element_density
while change >0.01 :
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
        filter_N[:,i]=NURBS
        density_basis[i]=Emin+(element_density[i]**penal)*(E0-Emin)
        K_E=density_basis[i]*K_E
        K_G=assemble(K_G,K_E,el_in,ncp)

    
    BC=BC_Switcher(CONTROL_POINTS,length,height,width,bc_disp)
    fixed_dof,load_dof,fixed_pts,load_pts=BC.indirect(option)
    reduced_k,reduced_F=apply_BC(K_G,F_E,fixed_dof,load_dof,load)

    U=np.matmul(np.linalg.inv(reduced_k),reduced_F)

    for j in fixed_dof:
        U=np.insert(U,j,0)
    
    F_E[load_dof]=load
    #print(U@F_E)
    
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
        element_density[k]=np.dot(node_density[k,:],NURBS)
        
        density_basis[k]=Emin+(element_density[k]**penal)*(E0-Emin)
        density_basis_dx[k]=-penal*(element_density[k]**(penal-1))*(E0-Emin)
        dof_index=np.sort(np.concatenate((el_in*dof,dof*el_in+1,dof*el_in+2)))
        compliance+=np.transpose(U[dof_index])@(density_basis[k]*K_E)@U[dof_index]
        dcompliance[k]=np.transpose(U[dof_index])@(density_basis_dx[k]*K_E)@U[dof_index]
    #compliance+=np.transpose(U)@(K_G)@U
    dv=np.ones(nel)
    #dcompliance=(1/np.maximum(1e-3,element_density)*DH)*np.matmul(H,(element_density*dcompliance))
    #delement_density=beta*np.exp(-beta*element_density)+np.exp(-beta)
    #dcompliance=np.matmul(H,(delement_density*dcompliance/DH))
    #dcompliance=np.matmul(H,(dcompliance/DH))
    #dv=np.matmul(H,dv/DH)
    #print(dcompliance)
    element_density=np.round(element_density,5)
    old_el_density=element_density
    #np.append(X,element_density)
    #element_density=np.matmul(H,element_density/DH)
    if loop>100:
        oc_disp=False
        bc_disp=False
        fil_disp=False
    if change<0.015:
        oc_disp=True
    constrain=nel*volume_frac
    if optimizer=='OC':
        dcompliance=(1/np.maximum(1e-3,element_density)*DH)*np.matmul(H,(element_density*dcompliance))

        element_density_updated= optimality_criteria(dcompliance,dv,element_density,constrain,H,DH,beta,oc_disp,g)
        element_density=element_density_updated
    if optimizer=='MMA':

        dcompliance=np.matmul(H,(dcompliance/DH))
        dv=np.matmul(H,dv/DH)

        v = (np.sum(element_density)/nel)-volume_frac
        #print(v)
        
        dv_dx  = (dv/ (volume_frac*nel))
        #print(dcompliance)
        #print(dv)
        dfun0=dcompliance
        dcon=dv_dx
        f0=compliance
        c0=v
        element_density_updated,Lower,Upper=Moving_asymptoes(dfun0,dcon,f0,c0,element_density,E1,E2,Lower,Upper,loop,nel,1,Xmin,Xmax,True)
        E2 = E1.copy()
        E1 = element_density.copy()
        element_density=np.round(element_density_updated[:,0],4)
    #Heavyside filter
    #filter_density=Heavyside_filter(element_density,beta,fil_disp)
    #print(element_density)
    #element_density=filter_density
    if (loop%50)==0:
        beta*=2
    CC.append(compliance)
    ii.append(loop)
    VV.append(np.mean(element_density))
    change=np.linalg.norm(element_density-old_el_density)
    #if optimizer=='MMA':
    #    element_density=np.matmul(H,element_density/DH)
    if loop==0:
        initial_compliance=compliance
    loop+=1
    print('     Iteration: %d       p :%f      optimizer:%s      OBJ=%f6       Vol=%f       CHG=%f6      '%(loop,penal,optimizer,compliance,np.mean(element_density),change)) 
    print('#'*width1) 
    print('='*width1)
    print('\n') 
    if loop >=max_iterations:
        print('~'*width1)
        fmt='{:^'+str(width1)+'}'
        print(fmt.format('Maximum iteration has been reached'))
        print('~'*width1)
        break
print('+'*width1) 

Mnd=np.dot(4*element_density,(1-element_density))/nel*100
print('                                        Measure of discreteness: ',Mnd)

element_density=np.round(element_density,1)
#element_density[element_density<0.1] = 0
#element_density[element_density>=volume_frac] = 1

print(np.mean(element_density))
print('\n                               Mass of the beam at 100% volume          :',(length*width*height)*density)
print('\n                               Optimized mass of the beam               :',(length*width*height)*density*np.mean(element_density))
energy_stored=np.dot(0.5,np.matmul(U.transpose(),np.matmul(U,K_G)))
print('\n                               Initial Compliance without optimization  :',initial_compliance)
final_compliance=np.dot(U,F_E)
print('\n                               Final Compliance with optimization       :',final_compliance,'\n')
'''
[1.   1.   1.   0.65 1.   1.   1.   0.7  1.   1.   1.   1.   1.   1.
 0.   0.   1.   1.   0.   0.   1.   1.   0.65 0.  ]
'''

print('%'*width1) 
print('\n Creating VTK file ')
from pyevtk.hl import gridToVTK

x=np.array(np.unique(CONTROL_POINTS[:,0]),dtype='float64')
y=np.array(np.unique(CONTROL_POINTS[:,1]),dtype='float64')
z=np.array(np.unique(CONTROL_POINTS[:,2]),dtype='float64')

volume = (element_density.reshape((nU, nV, nW),order='F'))

print(volume)
#comments = [ "comment 1", "comment 2" ]
gridToVTK("./cantilever_beam",x, y, z,cellData={"Volume":volume})

print('\n VTK file generated \n')
print('%'*width1)  


fig, ax = plt.subplots()
ax.plot(ii,CC)
ax.set_xlabel('iteration')
ax.set_ylabel('Compliance')
ax.set_xlim(0,ii[-1]+2)
#ax.set_ylim(0,CC)
ax.set_title('Compliance change in each iteration')
fig.savefig('ComplianceVSiteration.png')

fig, ax1 = plt.subplots()
ax1.plot(ii,VV)
ax1.set_xlabel('iteration') 
ax1.set_ylabel('Volume fraction')
ax1.set_xlim(0,ii[-1]+2)
ax.set_ylim(volume_frac-0.1,volume_frac+0.1)
ax1.set_title('Volume fraction change in each iteration')
fig.savefig('Volume_fractionVSiteration.png')
plt.show()


import pyvista as pv

pv.set_plot_theme("night")
mesh = pv.StructuredGrid()
# Set the coordinates from the numpy array
mesh.points = CONTROL_POINTS[:,:-2]
# set the dimensions
mesh.dimensions = [nx, ny, nz]
mesh.cell_arrays["volume"] =element_density
# and then inspect it!
'''
mesh.plot(show_edges=True,opacity='linear', cpos="xy")
'''


# read the data
# plot the data with an automatically created Plotter
#grid.plot(show_scalar_bar=True,opacity="linear",show_axes=True)


plotter = pv.Plotter(shape=(1, 2))
plotter.subplot(0, 0)
grid= mesh
#pv.set_plot_theme("night")
plotter.add_text("With opacity", font_size=10)
plotter.add_mesh(grid,opacity="linear")
plotter.show_bounds(all_edges=True)
#plotter.add_mesh(grid,opacity="linear",show_edges=True)
plotter.add_axes(interactive=True)

plotter.subplot(0, 1)
grid= mesh
#plotter.set_plot_theme("night")
plotter.add_text("Without opacity", font_size=10)
plotter.add_mesh(grid)
plotter.add_axes(interactive=True)
plotter.show_bounds(all_edges=True)
plotter.show(screenshot='canti.png')
