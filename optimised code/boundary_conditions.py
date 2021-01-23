import numpy as np
from Inputs import *
class BC_Switcher(object):

    def __init__(self,CONTROL_POINTS,length,height,width,bc_disp):
        self.CONTROL_POINTS=CONTROL_POINTS
        self.length=length
        self.width=width
        self.height=height
        self.bc_disp=bc_disp
        
    def indirect(self,i):
        method_name='number_'+str(i)
        method=getattr(self,method_name,lambda :'0')
        return method()
    
    
    def number_0(self):
        #nz should be odd
        if self.bc_disp:
            width1=120
            print('\n')
            print('='*width1)
            fmt='{:^'+str(width1)+'}'
            print(fmt.format(' Cantilever beam with load along the bottom edge of the free end'))
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==self.length and self.CONTROL_POINTS[index,1]==0] )
        load_indicies=np.sort((dof*load_nodes+1))


        return fixed_indicies,load_indicies,fixed_nodes,load_nodes
    def number_1(self):
        #nx and ny should be odd elements
        if self.bc_disp:
            width=120
            print('\n')
            print('='*width)
            fmt='{:^'+str(width)+'}'
            print(fmt.format(' Simple supported with load at bottom center'))
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if ((self.CONTROL_POINTS[index,0]==0 or self.CONTROL_POINTS[index,0]==self.length) and self.CONTROL_POINTS[index,1]==0 and (self.CONTROL_POINTS[index,2]==0 or self.CONTROL_POINTS[index,2]==self.width))])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if (self.CONTROL_POINTS[index,0]==self.length/2 and  self.CONTROL_POINTS[index,1]==0 and self.CONTROL_POINTS[index,2]==self.width/2 )] )
        load_indicies=np.sort((dof*load_nodes+1))


        return fixed_indicies,load_indicies,fixed_nodes,load_nodes

    def number_2(self):
        
        if self.bc_disp:
            width=120
            print('\n')
            print('='*width)
            fmt='{:^'+str(width)+'}'
            print(fmt.format(' Cantilever beam with load at end'))
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if (self.CONTROL_POINTS[index,0]==self.length and  (self.CONTROL_POINTS[index,1]==self.height or self.CONTROL_POINTS[index,1]==0)  and (self.CONTROL_POINTS[index,2]==self.width or self.CONTROL_POINTS[index,2]==0)) ] )
        load_indicies=np.sort((dof*load_nodes+1))


        return fixed_indicies,load_indicies,fixed_nodes,load_nodes

    def number_3(self):
        if self.bc_disp:
            width1=120
            print('\n')
            print('='*width1)
            fmt='{:^'+str(width1)+'}'
            print(fmt.format('Cantilever beam with point load at the free end (2d case loading at y=height and x=length) \n'))
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if (self.CONTROL_POINTS[index,0]==self.length and  self.CONTROL_POINTS[index,1]==0 and self.CONTROL_POINTS[index,2]==self.width/2 )] )
        load_indicies=np.sort((dof*load_nodes+1))
        return fixed_indicies,load_indicies,fixed_nodes,load_nodes
    
    def number_4(self):
        if self.bc_disp:
            width=120
            print('\n')
            print('='*width)
            fmt='{:^'+str(width)+'}'
            print(fmt.format('Cantilever beam with point load at the free end'))
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if (self.CONTROL_POINTS[index,0]==self.length and  self.CONTROL_POINTS[index,1]==self.height and  self.CONTROL_POINTS[index,2]==0)])
        load_indicies=np.sort((dof*load_nodes+1))


        return fixed_indicies,load_indicies,fixed_nodes,load_nodes



'''
s=Switcher()
s.indirect(2) 
print(CONTROL_POINTS)
print(s.indirect(0) )

'''
'''
length=48
height=5
width=10
'''
#provide the number of elements in each direction

'''
nx=21
ny=8
nz=11
'''
#Provide the  information for knot vector in each direction
'''

N=nx+1
P=ny+1
Q=nz+1

XI_DEGREE=1
ETA_DEGREE=1
NETA_DEGREE=1

C=Inputs(length,height,width,nx,ny,nz,XI_DEGREE,ETA_DEGREE,NETA_DEGREE)

CONTROL_POINTS=C.crtpts_coordinates()
print(CONTROL_POINTS)
bc_disp=False
BC=BC_Switcher(CONTROL_POINTS,length,height,width,bc_disp)
print(BC.indirect(1))
'''