import pytest
from optimization import optimality_criteria,Moving_asymptoes

from numpy import array,round,array_equiv,linalg,equal,diag,transpose,ones,zeros,dot,size,maximum,concatenate


def test__optimality_criteria(): 
    '''
    finding optimality conditions based on KKT conditions
    '''
    def fun(x):
        '''
        linear function
        f=2*x1+x2-2=0
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
    while iteration<50:
        iteration+=1
        xold=x
        df=dfun(x)
        dc=dcon(x)
        constrain=1
        xlist=optimality_criteria(df,dc,x,constrain)
        x=xlist
        if linalg.norm(x-xold)<1e-4:
            exit
    x_exact=array([0.5,0.5])
    x=round(x,2)
    value=fun(x)
    #print(value)
    assert (array_equiv(x,x_exact) and equal(value,0.0)) is True
    

def test__Moving_Asymptoes():

    def quad_function(x):
        return x[0]**2+x[1]**2+x[2]**2

    def dfunction(x):
        df_dx=array([2*x[0],2*x[1],2*x[2]])
        return df_dx

    def constrain(x):
        c1=(x[0]-5)**2+(x[1]-2)**2+(x[2]-1)**2-9
        c2=(x[0]-3)**2+(x[1]-4)**2+(x[2]-3)**2-9
        C=array([c1,c2])
        return C

    def dconstrains(x):
        c11=2*(x[0]-5)
        c12=2*(x[1]-2)
        c13=2*(x[2]-1)
        c21=2*(x[0]-3)
        c22=2*(x[1]-4)
        c23=2*(x[2]-3)
        dc=array([[c11,c12,c13],[c21,c22,c23]])
        return dc
    
    xval = array([[4,3,2]]).T
    nx=len(xval)

    xmin = zeros(nx)
    xmax = 5*ones(nx)
    L = xmin
    U = xmax

    loop=0
    max_loop=50

    while loop<=max_loop:
        vel=2
        nel=3
        xold=xval[:,loop]
        f0=quad_function(xold)
        df0=dfunction(xold)
        c0=constrain(xold)
        dc0=dconstrains(xold)
        x,L,U=Moving_asymptoes(df0,dc0,f0,c0,xval,L,U,loop,nel,vel,xmin,xmax,True)
        xval=concatenate((xval,x),axis=1)
        #print(xval)
        loop+=1
        xnew=xval[:,loop]

        #print(xold-xnew)
        if abs(linalg.norm(xold-xnew))<1e-4:
            f0=quad_function(xnew)
            df0=dfunction(xnew)
            c0=constrain(xnew)
            dc0=dconstrains(xnew)
            print(x,df0,dc0,f0,c0)
            break
    x=round(xnew,2)
    x_exact=round(array([2.0175, 1.7800, 1.237]),2)
    assert array_equiv(x,x_exact) is True

    #assert (0==0) is True

