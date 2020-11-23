import numpy as np
import pytest


class Inputs():
    '''
    A geometry class which
    '''
    def __init__(self,length=2,height=2,width=2,nx=2,ny=2,nz=2,xidegree=2,etadegree=2,netadegree=2):
        self.length=length
        self.height=height
        self.width=width
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.n=(nx+1)
        self.p=(ny+1)
        self.q=(nz+1)
        self.xidegree=xidegree
        self.etadegree=etadegree
        self.netadegree=netadegree



    def crtpts_coordinates(self):
        '''
        A function which return the control point co-ordinates
        '''
        beam_coordinates=[]
        index=0
        for i in range(self.nz+1):
            for j in range(self.ny+1):
                for k in range(self.nz+1):
                    beam_coordinates.append([(k)*(self.length/self.nx),j*(self.height/self.ny),i*(self.width/self.nz),1])
                    index+=1
        return np.array(beam_coordinates)
    


    def knot_vector(self,n,degree):
        '''
        This function return knotvector based on the number of control points and degree of Bspline
        Parameters
        ----------
        control_points : int
            number of points along with weights.
        degree : int
            order of B-splines  0-constant, 1-linear, 2-quadratic, 3-cubic.

        Returns
        -------
        KNOTVECTOR - an array containing knots based on control points
        '''
        size=n+degree+1
        self.knotvector=np.zeros(size)
        if n >=(degree):
            for i in range(size):
                if i<=degree:
                    self.knotvector[i]=0
                elif i>degree and i<(size-degree):
                    value=i-degree
                    self.knotvector[i]=value
                else:
                    self.knotvector[i]=n-degree
        elif n==2 and degree==1:
            self.knotvector=[0,0,1,1]
        elif n==2 and degree==2:
            self.knotvector=[0,0,0,1,1,1]
        else:
            raise Exception('Require more control points to define the knotvector')
        self.knotvector=self.knotvector/self.knotvector[-1]
        return self.knotvector


    def xi_knotvector(self):
        self.xikntvtr=self.knot_vector(self.n,self.xidegree)
        return self.xikntvtr

    def eta_knotvector(self):
        self.etakntvtr=self.knot_vector(self.p,self.etadegree)
        return self.etakntvtr
    
    def neta_knotvector(self):
        self.netakntvtr=self.knot_vector(self.q,self.netadegree)
        return self.netakntvtr


