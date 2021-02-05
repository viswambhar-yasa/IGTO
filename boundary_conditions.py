import numpy as np
from test_preprocessing import Inputs
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
            TBLACK =  '\033[30;1;47m' # Green Text
            ENDC = '\033[m'
            print(TBLACK+fmt.format(' Cantilever beam with load along the bottom edge of the free end')+ENDC)
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
            TBLACK =  '\033[30;1;47m' # Green Text
            ENDC = '\033[m'
            fmt='{:^'+str(width)+'}'
            print(TBLACK+fmt.format(' Simple supported with load at bottom center')+ENDC)
        dof=3
        fixed_nodes1=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if ((self.CONTROL_POINTS[index,0]==0 ) and self.CONTROL_POINTS[index,1]==0 and (self.CONTROL_POINTS[index,2]==0 or self.CONTROL_POINTS[index,2]==self.width))])
        fixed_indicies1=np.sort(np.concatenate((fixed_nodes1*dof,dof*fixed_nodes1+1,fixed_nodes1*dof+2)))
        fixed_nodes2=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if ((self.CONTROL_POINTS[index,0]==self.length ) and self.CONTROL_POINTS[index,1]==0 and (self.CONTROL_POINTS[index,2]==0 or self.CONTROL_POINTS[index,2]==self.width))])
        fixed_indicies2= np.sort(np.concatenate((dof*fixed_nodes2+1,fixed_nodes2*dof+2)))
        fixed_nodes=np.sort(np.concatenate((fixed_nodes1,fixed_nodes2)))
        fixed_indicies=np.sort(np.concatenate((fixed_indicies1,fixed_indicies2)))
        load_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if (self.CONTROL_POINTS[index,0]==self.length/2 and  self.CONTROL_POINTS[index,1]==self.height and self.CONTROL_POINTS[index,2]==self.width/2 )] )
        load_indicies=np.sort((dof*load_nodes+1))


        return fixed_indicies,load_indicies,fixed_nodes,load_nodes

    def number_2(self):
        
        if self.bc_disp:
            width=120
            print('\n')
            print('='*width)
            TBLACK =  '\033[30;1;47m' # Green Text
            ENDC = '\033[m'
            fmt='{:^'+str(width)+'}'
            print(TBLACK+fmt.format(' Cantilever beam with load at end')+ENDC)
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
            TBLACK =  '\033[30;1;47m' # Green Text
            ENDC = '\033[m'
            fmt='{:^'+str(width1)+'}'
            print(TBLACK+fmt.format('Cantilever beam with point load at the free end (2d case loading at y=height and x=length) ')+ENDC)
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
            TBLACK =  '\033[30;1;47m' # Green Text
            ENDC = '\033[m'
            print(TBLACK+fmt.format('Cantilever beam with two forces at either end')+ENDC)
        dof=3
        fixed_nodes=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if self.CONTROL_POINTS[index,0]==0])
        fixed_indicies=np.sort(np.concatenate((fixed_nodes*dof,dof*fixed_nodes+1,fixed_nodes*dof+2)))
        load_nodes=[]
        load_indicies=[]
        load_nodes1=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if (self.CONTROL_POINTS[index,0]==self.length and  self.CONTROL_POINTS[index,1]==0 and  self.CONTROL_POINTS[index,2]==self.width/2)])
        load_nodes.append(load_nodes1)
        load_indicies.append(np.sort((dof*load_nodes1+1)))
        load_nodes1=np.array([index for index,j in enumerate(self.CONTROL_POINTS) if (self.CONTROL_POINTS[index,0]==self.length and  self.CONTROL_POINTS[index,1]==self.height and  self.CONTROL_POINTS[index,2]==self.width/2)])
        load_nodes.append(load_nodes1)
        load_indicies.append(np.sort((dof*load_nodes1+1)))
        return fixed_indicies,load_indicies,fixed_nodes,load_nodes


