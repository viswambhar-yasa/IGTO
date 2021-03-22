#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#------------------------------------#
#A python test file to check optimizers 
# command to run all test cases
# pytest test_optimization.py
# -----------------------------------# 
import pytest
from optimization import optimality_criteria,Moving_asymptoes
from numpy import array,round,array_equiv,linalg,equal,zeros,ones,all


def test__optimality_criteria_simple_function(): 
    '''

    
    finding optimality conditions using penality method for constrains like the 
     sum(x0,x1,x2,x4,...)=1
     and 0<=x<=1
    '''
    def fun(x):
        '''
        linear function
        f(x1,x2)=x1+x2-2=0
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
    residual=0
    while iteration<50:
        iteration+=1
        xold=x
        df=dfun(x)
        dc=dcon(x)
        constrain=1
        xlist,residual=optimality_criteria(df,dc,x,constrain)
        x=xlist
        xplot[iteration,:]=x
        if linalg.norm(x-xold)<1e-4:
            exit
    x_exact=array([0.5,0.5])
    error=all(abs(x-x_exact)<1e-2)
    o=0
    if error:
        o=1
    assert (o==1) is True

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
    residual=0
    while iteration<50:
        iteration+=1
        xold=x
        df=dfun(x)
        dc=dcon(x)
        constrain=1
        xlist,residual=optimality_criteria(df,dc,x,constrain)
        x=xlist

        xplot[iteration,:]=x
        #print(xplot)
        if linalg.norm(x-xold)<1e-3:
            exit
    x_exact=array([0,1])
    x=round(x,3)
    error=all(abs(x-x_exact)<1e-2)
    o=0
    if error:
        o=1
    assert (o==1) is True  



def test__MMA_literature_equation():
    
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
    x=array([4,3,2])
    f0=quad_function(x)
    df0=dfunction(x)
    c0=constrain(x)
    dc0=dconstrains(x)
    #print(x,df0,dc0,f0,c0)

    xval = array([4,3,2])
    print(xval)
    nx=len(xval)

    xmin = zeros(nx)
    xmax = 5*ones(nx)

    L = xmin
    U = xmax

    loop=0
    max_loop=50
    x1=xval
    x2=xval
    l=0
    all_values=[]
    while loop<=max_loop:
        vel=2
        nel=3
        f0=quad_function(xval)
        df0=dfunction(xval)
        c0=constrain(xval)
        dc0=dconstrains(xval)
        x,L,U=Moving_asymptoes(df0,dc0,f0,c0,xval,x1,x2,L,U,loop,nel,vel,xmin,xmax,False)
        print(x.T)

        x2 = x1.copy()
        x1 = xval.copy()
        x3=x.T
        xval=x3[0]
        print(xval)
        all_values.append(x3)
        #print(xold-xnew)

        l+=1
        if abs(linalg.norm(x1-xval))<0.0001:
            f0=quad_function(xval)
            df0=dfunction(xval)
            c0=constrain(xval)
            dc0=dconstrains(xval)
            print(x3,df0,dc0,f0,c0)
            break
    x_output=x3[0]
    print(x_output)
    x_exact=array([2.0175, 1.7800, 1.237])
    error= all(abs(x_exact-x_output)<1e-3)
    o=0
    if error:
        o=1
    assert (o==1) is True
