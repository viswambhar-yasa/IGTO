import pytest
from optimization import optimality_criteria
from numpy import array,round,array_equiv,linalg,equal


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
        xlist=optimality_criteria(df,dc,x,constrain,oc_disp=False)
        x=xlist
        if linalg.norm(x-xold)<1e-4:
            exit
    x_exact=array([0.5,0.5])
    x=round(x,2)
    value=fun(x)
    assert (array_equiv(x,x_exact) and equal(value,0.0)) is True
    