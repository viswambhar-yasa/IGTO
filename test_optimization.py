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
    FUNCTIONAL TESTING
    AIM : To check the performance of optimal criteria function by solving a simple equation with simple constrains
            Linear function
            f(x1,x2)=x1+x2-2=0
            subjected to 
            inequality constrained
            constrain=sum(x1+x2)<=1
            x1<1
            x2<1
    Expected result : [0.5,0.5]
    Remarks :passed sucessfully
 
    '''
    def fun(x):
        '''
        linear function
        f(x1,x2)=x1+x2-2=0
        '''
        f=x[0]+x[1]-1
        return f

    def dfun(x):
        '''
        derivatives of function
        df/dx1=1
        fd/dx2=1
        '''
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
    #initial guess to obtain optimal solution
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
        #updated x variables obtained from optimal criteria (penality based method)
        xlist,residual=optimality_criteria(df,dc,x,constrain)
        x=xlist
        xplot[iteration,:]=x
        #exit conditon
        if linalg.norm(x-xold)<1e-4:
            exit
    x_exact=array([0.5,0.5])
    error=all(abs(x-x_exact)<1e-2)
    o=0
    if error:
        o=1
    assert (o==1) is True

def test__optimality_criteria_quadratic_function():
    '''
    FUNCTIONAL TESTING
    AIM : To check the performance of optimal criteria function by solving a quadratic equation with simple constrains
            Linear function
            f=(x0+x1)^2-2*(x0-x1)^2-1
            subjected to 
            inequality constrained
            constrain=sum(x1+x2)<=1
            x1<1
            x2<1
    Expected result : [0.5,0.5]
    Remarks :passed sucessfully
    '''
    def fun(x):
        '''
            quadratic function
            f=(x0+x1)^2-2*(x0-x1)^2-1
        '''
        f=(x[0]+x[1])**2-2*(x[0]-x[1])**2-1
        return f

    def dfun(x):
        '''
        derivatives of function
        df/fx0=-2*x[0]+6*x[1]
        df/dx1=6*x[1]-2*x[1]
        '''
        #df_dx0=2*(x[1])+2
        #df_dx1=4*x[1]+2*x[0]
        df_dx0=-2*x[0]+6*x[1]
        df_dx1=-2*x[1]+6*x[0]
        df=array([df_dx0,df_dx1])
        return df

    def dcon(x):
        '''
        derivatives of constrains
        '''
        dc_dx0=-1
        dc_dx1=-1
        dc=array([dc_dx0,dc_dx1])
        return dc
    #initial guess to obtain optimal solution
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
        #updated x variables obtained from optimal criteria (penality based method)
        xlist,residual=optimality_criteria(df,dc,x,constrain)
        x=xlist

        xplot[iteration,:]=x
        #print(xplot)
        #exit conditon
        if linalg.norm(x-xold)<1e-3:
            exit
    x_exact=array([0.5,0.5])
    x=round(x,3)
    print(x)
    error=all(abs(x-x_exact)<1e-2)
    o=0
    if error:
        o=1
    assert (o==1) is True  



def test__MMA_literature_equation():
    '''
    FUNCTIONAL TESTING
    AIM : To check the performance of method of moving asymptotes function by solving a quadratic equation with complex constrains.
        The equation and constrains and expected outputs are obtained from MMA paper.
            A quadratic function
            f(x0,x1,x2)=x0^2+x1^2+x3^2
            subjected to 
                c1=(x0-5)**2+(x1-2)**2+(x2-1)**2-9
                c2=(x0-3)**2+(x1-4)**2+(x2-3)**2-9
                0 =< x(j) =< 5; for j = 1; 2; 3
    Expected result : [2.0175, 1.7800, 1.237]
    Remarks :passed sucessfully
    '''
    def quad_function(x):
        '''
        A quadratic function
        f(x0,x1,x2)=x0^2+x1^2+x3^2
        '''
        return x[0]**2+x[1]**2+x[2]**2

    def dfunction(x):
        '''
        derivative of the function
        df_dx0=2*x0
        df_dx1=2*x1
        df_dx2=2*x2
        '''
        df_dx=array([2*x[0],2*x[1],2*x[2]])
        return df_dx

    def constrain(x):
        '''
        Constrains of the function
        '''
        c1=(x[0]-5)**2+(x[1]-2)**2+(x[2]-1)**2-9
        c2=(x[0]-3)**2+(x[1]-4)**2+(x[2]-3)**2-9
        C=array([c1,c2])
        return C

    def dconstrains(x):
        '''
        derivative of constrains w.r.t variables x0,x1,x2
        '''
        c11=2*(x[0]-5)
        c12=2*(x[1]-2)
        c13=2*(x[2]-1)
        c21=2*(x[0]-3)
        c22=2*(x[1]-4)
        c23=2*(x[2]-3)
        dc=array([[c11,c12,c13],[c21,c22,c23]])
        return dc
    #initial guess    
    x=array([4,3,2])
    f0=quad_function(x)
    df0=dfunction(x)
    c0=constrain(x)
    dc0=dconstrains(x)
    #print(x,df0,dc0,f0,c0)

    xval = array([4,3,2])
    print(xval)
    nx=len(xval)
    #initializing the parameter required for MMA 
    #maximum and minimum of the variables
    xmin = zeros(nx)
    xmax = 5*ones(nx)

    L = xmin
    U = xmax

    loop=0
    max_loop=50
    #MMA requires previous two iteration values
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
        #the optimal solution is obtained from moving asymptoes (gradient based method)
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
        if abs(linalg.norm(x1-xval))<0.0001: #exit condition
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
