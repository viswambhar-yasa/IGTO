import pytest
from optimization import optimality_criteria
from numpy import array,round,array_equiv,linalg,equal,zeros


def test__optimality_criteria_simple_function(): 
    '''
    finding optimality conditions using penality method for constrains like the 
     sum(x0,x1,x2,x4,...)=1
     and 0<=x<=1
    '''
    def fun(x):
        '''
        linear function
        f=x1+x2-2=0
        '''
        f=x[0]+x[1]-1
        return f

    def dfun(x):
        df_dx1=1
        df_dx2=1
        df=array([df_dx1,df_dx2])
        return df

    def dcon(x):
        '''
        inequality constrained
        constrain=sum(x1+x2)<=1
        x1<1
        x2<1
        '''
        dc_dx1=-1
        dc_dx2=-1
        dc=array([dc_dx1,dc_dx2])
        return dc

    initial_x=array([1,1])

    iteration=1
    x=initial_x
    xplot=zeros((51,len(initial_x)))
    while iteration<50:
        iteration+=1
        xold=x
        df=dfun(x)
        dc=dcon(x)
        constrain=1
        xlist=optimality_criteria(df,dc,x,constrain)
        x=xlist
        xplot[iteration,:]=x
        if linalg.norm(x-xold)<1e-4:
            exit
    x_exact=array([0.5,0.5])
    assert (all(abs(x-x_exact)<1e-2)) is True

def test__optimality_criteria_quadratic_function():
    
    def fun(x):
            '''
            quadratic function
            f=(x0+x1)^2-2*(x0-x1)
            '''
            f=(x[0]+x[1])**2-2*(x[0]-x[1])-1
            return f

    def dfun(x):
        df_dx0=2*(x[1])
        df_dx1=4*x[1]+2*x[0]
        df=array([df_dx0,df_dx1])
        return df

    def dcon(x):
        dc_dx0=-1
        dc_dx1=-1
        dc=array([dc_dx0,dc_dx1])
        return dc

    initial_x=array([1,1])
    iteration=1
    x=initial_x
    xplot=zeros((51,len(initial_x)))
    while iteration<50:
        iteration+=1
        xold=x
        df=dfun(x)
        dc=dcon(x)
        constrain=1
        xlist=optimality_criteria(df,dc,x,constrain)
        x=xlist

        xplot[iteration,:]=x
        print(xplot)
        if linalg.norm(x-xold)<1e-3:
            exit
    x_exact=array([0,1])
    x=round(x,3)
    assert all(abs(x-x_exact)<1e-2) is True  


