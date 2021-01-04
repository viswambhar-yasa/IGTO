import numpy as np
from Inputs import *
class BC_Switcher(object):

    def __init__(self,CONTROL_POINTS,length,height,width):
        self.CONTROL_POINTS=CONTROL_POINTS
        self.length=length
        self.width=width
        self.height=height
        
    def indirect(self,i):
        method_name='number_'+str(i)
        method=getattr(self,method_name,lambda :'0')
        return method()
    
    
    def number_0(self):

        print('\n Cantilever beam with point load at the free end \n')
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==self.length and self.CONTROL_POINTS[index,1]==self.height and self.CONTROL_POINTS[index,2]==self.width] )
        load_indicies=np.sort((dof*load_nodes+1))


        return fixed_indicies,load_indicies,fixed_nodes,load_nodes
    def number_1(self):

        print('\n Cantilever beam with pressure at top \n')
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,1]==self.height and self.CONTROL_POINTS[index,0]!=0 ] )
        load_indicies=np.sort((dof*load_nodes+1))


        return fixed_indicies,load_indicies,fixed_nodes,load_nodes

    def number_2(self):
        print('\n Cantilever beam with pressure at top \n')
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==self.length and self.CONTROL_POINTS[index,1]==0 and self.CONTROL_POINTS[index,2]==0  ] )
        load_indicies=np.sort((dof*load_nodes+1))


        return fixed_indicies,load_indicies,fixed_nodes,load_nodes

    def number_3(self):
        print('\n 2d case loading at y=height and x=length \n')
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if (self.CONTROL_POINTS[index,0]==self.length and  self.CONTROL_POINTS[index,1]==self.height and self.CONTROL_POINTS[index,2]==0 )] )
        load_indicies=np.sort((dof*load_nodes+1))
        return fixed_indicies,load_indicies,fixed_nodes,load_nodes
    
    def number_4(self):

        print('\n Cantilever beam with point load at the free end \n')
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))

        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if (self.CONTROL_POINTS[index,0]==self.length and  self.CONTROL_POINTS[index,1]==self.width and  self.CONTROL_POINTS[index,2]==float(self.height/2))] )
        load_indicies=np.sort((dof*load_nodes+1))


        return fixed_indicies,load_indicies,fixed_nodes,load_nodes



'''
s=Switcher()
s.indirect(2) 
print(CONTROL_POINTS)
print(s.indirect(0) )

'''

length=1
height=0.1
width=0.1

'''
provide the number of elements in each direction

'''
nx=6
ny=6
nz=7

'''
Provide the  information for knot vector in each direction
'''

N=nx+1
P=ny+1
Q=nz+1

XI_DEGREE=1
ETA_DEGREE=1
NETA_DEGREE=1

C=Inputs(length,height,width,nx,ny,nz,XI_DEGREE,ETA_DEGREE,NETA_DEGREE)
'''
CONTROL_POINTS=C.crtpts_coordinates()
print(CONTROL_POINTS)
BC=BC_Switcher(CONTROL_POINTS,length,height,width)
print(BC.indirect(4))
'''