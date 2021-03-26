#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#---------------------------------------------------------------------------------------#
#A python file where Iso-geometric analysis is performed along with Toplogy optimization
# --------------------------------------------------------------------------------------# 

import numpy as np
import pytest

def Knearestneighbours(rmin,nelx,nely,nelz):
    '''

    The code has been taken from An efficient 3D topology optimization code written in Matlab section
    To perform sensitivity analysis, we require weight function.This functions provides the weight factor 


    Parameters
    ----------
    rmin : int
        The minimum radius in which the elements are influenced.
    nelx,nely,nelz : int
        No of element in x ,y and z direction.

    Returns
    -------
    H : array
        A 2d array containing the weight factor of each element.
    DH : array
        Sum of the weight factor (change in weight factors).

    '''
    nel=nelx*nely*nelz #total number of elements
    #initializing the dimensions of weight factor
    H=np.zeros((nel,nel))
    #looped over element along x y and z 
    for i in range(nelz):
        for j in range(nely):
            for k in range(nelx):
                r=int(i*nelx*nely+j*nelx+k)
                imin=int(max(i-(np.ceil(rmin)-1),0))
                imax=int(min(i+np.ceil(rmin),nelz))
                jmin=int(max(j-(np.ceil(rmin)-1),0))
                jmax=int(min(j+np.ceil(rmin),nely))
                kmin=int(max(k-(np.ceil(rmin)-1),0))
                kmax=int(min(k+np.ceil(rmin),nelx))
                for ii in range(imin,imax): 
                    for jj in range(jmin,jmax):
                        for kk in range(kmin,kmax):
                            c= int(ii*nelx*nely+jj*nelx+kk)
                            distance=np.sqrt((i-ii)**2+(j-jj)**2+(k-kk)**2)
                            H[r,c]=np.maximum(0,(rmin-distance))        # based on equ.4.33
    H=np.array(H)
    DH=np.sum(H,1)
    return H,DH  


def optimality_criteria(dfun,dCon,initial_X,constrain,residual=None,H=None,DH=None,oc_disp=True,g=1,tol=0.01,move=0.1,neta=0.5,initial_value=0,final_value=1e09):
    '''
    Optimality criteria is a penality based method to solve constrain based probelem.


    Parameters
    ----------
    dfun : array
        derivatives of the function which are to be optimized.
    dCon : array
        derivatives of the constrains .
    initial_X : array
        initial set of values which are required to solve the problem.
    constrain : array
                set of constrains.
    H : array
        weight factor required for sensitivity analysis.
    DH : array
        change in weight factor.
    oc_disp : bool, optional
        If TRUE the values at each iteration is printed. The default is True.
    g : float, optional
        penality value to remove gray scale values. The default is 1.
    tol : float, optional
        Tolerance . The default is 0.01.
    move : float, optional
        the step which has to be taken while find the optimization. The default is 0.1.
    neta : float, optional
        The power of the constrained problem. The default is 0.5.
    initial_value : int, optional
        The initial range to find lagarian multiplier. The default is 0.
    final_value : TYPE, optional
        The final range to find lagarian multiplier. The default is 1e09.

    Returns
    -------
    X_new : array 
        The optimized values obtained after optimization.

    Test case
    ---------
    test command -
         

    '''
    #initializing the values of optimized value
    X_new=np.zeros(len(initial_X))
    X=np.array(initial_X)
    i=0
    max_iteration=150
    while i<=max_iteration:
        convergence_criteria=(final_value-initial_value)/(final_value+initial_value)

        #lagarian multiplier
        lamdba=0.5*(initial_value+final_value)
        #Bey
        Be=((-dfun/(dCon*lamdba))**neta)**g #based on equ. 4.38
        X_new[:]= np.maximum(0.0,np.maximum(X-move,np.minimum(1.0,np.minimum(X+move,X*Be)))) #based on equ. 4.37
        if H is not None:
            X_new=np.matmul(H,X_new/DH) #based on equ.4.32
        #  bi-section algorthim updating lagaranian multiplier
        #multiple exit conditions to improve the functionality of the OC
        if residual is None:
                update_condition=np.sum(X_new)
                if update_condition>=constrain:
                    initial_value=lamdba    
                else:
                    final_value=lamdba
        else:
            update_condition=residual+np.sum((dCon*(X_new-X)))
            if update_condition>=0:
                initial_value=lamdba    
            else:
                final_value=lamdba
        
        
        i=i+1
        TBLUE = '\033[34;1m'
        TRED='\033[31;1m'
        ENDC = '\033[m'
        #exit condition
        if convergence_criteria<=tol :
            width=120
            if oc_disp:
                print('-'*width)
                fmt='{:^'+str(width)+'}'
                print(TBLUE+fmt.format('Optimiality Criterian')+ENDC)
                print('.'*width)
                print(TBLUE+'     exit_iter: %d    gray_filter(g):%f    lamdba:%f12      tot_vol:%f4     max_volume:%f4' %(i,g,lamdba,update_condition,constrain)+ENDC)  
                print('.'*width)
            return X_new,update_condition
        if i==(max_iteration-1):
            width=120
            if oc_disp:
                print(TRED+'maximum iteration has reached, didnot converger'+ENDC)
                print('.'*width)
                fmt='{:^'+str(width)+'}'
                print(TBLUE+fmt.format('Optimiality Criterian')+ENDC)
                print('.'*width)
                print(TBLUE+'     exit_iter: %d        lamdba:%f12        tot_vol:%f4       max_volume:%f4' %(i,lamdba,update_condition,constrain)+ENDC)  
                print('*'*width)
            return X_new,update_condition



def Moving_asymptoes(dfun0,dcon,f0,c0,x0,x1,x2,L,U,loop,nel,vel,Xmin,Xmax,ma_disp,m=0.2,m_tol=0.1):
    '''
    Based on the MMA paper by Krister Svanberg whose formulation are used to build this function
    MMA and GMMA by Krister Svanberg
    Optimization and Systems Theory,
    KTH, Stockholm, Sweden.
    The equation are provided in the paper and are referred in this function accordingly.
    This is the main function which performs active set strategy based constrained problem optimization (KKT condition are satisfied using Newton raphson)

    Parameters
    ----------
    dfun0 : array
        Contains the derivatives of the function.
    dcon : array
        Contains the derivatives of the constrains.
    f0 : array
        inital value of the function for value x0.
    c0 : array
        initial values of the constraind for values x0.
    x0 : array
        Initial guess values and later the current iteration values of X.
    x1 : array
        X values from previous iteration (k-1).
    x2 : array
        X values from previous iteration (k-2).
    L : array
        Lower limit values.
    U : array
        Upper limit values.
    loop : int
        Current iteration number.
    nel : int 
        Number of variables in the function. [x1,x2,x3,...]
    vel : int
        Number of constrains.
    Xmin : array
        The minimum values of  X.
    Xmax : array
        The maximum value of X.
    ma_disp : bool
        If TRUE , prints the values.
    m : float, optional
        The length of the step movement. The default is 0.2.
    m_tol : bool, optional
        Tolerance for MMA. The default is 0.1.

    Returns
    -------
    x  : array
            Optimized values of the array which satisfys KKT condition
    Lower : array
                The lower limit for the values of X
    Upper : array
            Upper limit for  the values of X

    Test case
    ---------
    Test command -
                pytest test_optimization.py::test__MMA_literature_equation

    '''
    
    
    def Asymptoes(loop=loop,x0=x0,x1=x1,x2=x2,n=nel,Xmin=Xmin,Xmax=Xmax,L=L,U=U,move_tol_low=0.01,move_tol_high=10):
        '''
        The lower and upper bound are calculated based on the X values in previous 2 iteration.

        Parameters
        ----------
        loop : int
            The current number of the iteration.
        x0,x1,x2 :array
                The values at loop,loop-1 and loop-2 iterations
        n : int
            Number of variable in x 
        Xmin,Xmax : array
                    The minimum and maximum values of X 
        L,U : array
            The lower and upper bound of x obtained from previous iteration
        move_tol_low : float, optional
            tolerance for lower limit . The default is 0.01.
        move_tol_high : int, optional
            tolerance for upper limit . The default is 10.

        Returns
        -------
        lower_value : array
            The Lower limit for the values of X[x1,x2,x3,...].
        Upper : array
            The Upper limit for the values of X[x1,x2,x3,...].
        x0 : array
            The initial values of X[x1,x2,x3,.....].

        '''

        lower_value=np.ones(n)
        upper_value=np.ones(n)
        #iteration less then 3
        if loop<=2:
            lower_value=x0-0.5*(Xmax-Xmin)  #based on equ. 3.11 MMA paper
            upper_value=x0+0.5*(Xmax-Xmin)  #based on equ. 3.11 MMA paper
            #print(Lower,Upper,x)
            return lower_value,upper_value,x0
        else: #iteration greater then 3
            xval=x0  #current iteration value of x
            xold1=x1 #pervious iteration value of x (k-1)
            xold2=x2 #pervious iteration value of x (k-2)
            zzz = (xval-xold1)*(xold1-xold2)  #based on equ. 3.13 MMA paper
            gamma = np.ones(nel)
            gamma[np.where(zzz>0)] = 1.2      #based on equ. 3.13 MMA paper
            gamma[np.where(zzz<0)] = 0.7      #based on equ. 3.13 MMA paper
            lower_value = xval-gamma*(xold1-L) #based on equ. 3.12 MMA paper
            #print(Lower)
            upper_value = xval+gamma*(U-xold1) #based on equ. 3.12 MMA paper
            #print(Upper)
            lmin = xval-move_tol_high*(Xmax-Xmin) #based on equ. 3.14 MMA paper
            lmax = xval-move_tol_low*(Xmax-Xmin)  #based on equ. 3.14 MMA paper
            umin = xval+move_tol_low*(Xmax-Xmin)  #based on equ. 3.14 MMA paper
            umax = xval+move_tol_high*(Xmax-Xmin) #based on equ. 3.14 MMA paper
            lower_value = np.maximum(lower_value,lmin)  #based on equ. 3.6 MMA paper
            lower_value = np.minimum(lower_value,lmax)  #based on equ. 3.6 MMA paper
            upper_value = np.minimum(upper_value,umax)  #based on equ. 3.7 MMA paper
            upper_value = np.maximum(upper_value,umin)  #based on equ. 3.7 MMA paper
            return lower_value,upper_value,x0
            ##smaple output as : [ 1.5  0.5 -0.5], [6.5 5.5 4.5], [4 3 2]
        
    
    def objective_constrains(x,Lower,Upper,dfun0=dfun0,f0=f0,Xmin=Xmin,Xmax=Xmax,tol=1e-5,tol1=0.001,tol2=1.001):
        '''
        Approximation of the function is built using function and it's derivatives p and q 
        Based on equ. 4.2 MMA paper

        Parameters
        ----------
        x : array
            The intial value of variables.
        Lower,Upper : array
            The lower and upper limit of the variables.
        fo,dfun0 : array, optional
                the function value and derivative of the function for initial values of variable
        Xmin,Xmax : array, optional
                        The maximum and minimum od the variables X
        tol : float, optional
            Tolerance for objective constrain. The default is 1e-5.
        tol1,tol2 : float, optional
                        Tolerance for objective constrain . The default is 0.001 and 1.001.

        Returns
        -------
        p : array
            function approximation equ. 4.2 in MMA paper.
        q : array
            function approximation equ. 4.2 in MMA paper
        '''
        df_postive=np.maximum(dfun0,0)   #Based on 4.3 in MMA paper # positive values of derivative function
        df_negative=np.maximum(-dfun0,0) #Based on 4.3 in MMA paper #negative values of derivative function
        #print(df_postive)
        Xd_inv=1/Xmax-Xmin  
        UX=Upper-x
        
        LX=x-Lower
        #
        p0=tol2*df_postive+tol1*df_negative+(tol*Xd_inv) 
        #print('p0',p0)
        p=np.array([(UX**2)*p0]).T  #Based on equ, 4.3 in MMA paper
        q0=tol1*df_postive+tol2*df_negative+(tol*Xd_inv) 
        q=np.array([(LX**2)*q0]).T  #Based on 4.4 in MMA paper

        #print('p',p)
        #print('q',q)    
        return p,q

    def Minimizer_constrains(x,Lower,Upper,dcon=dcon,m=vel,c0=c0,Xmin=Xmin,Xmax=Xmax,tol=1e-5,tol1=0.001,tol2=1.001):
        '''
        Approximation of the function is built using constrains and it's derivatives P and Q
        Based on equ. 4.2 MMA paper

        Parameters
        ----------
        x : array   
            Optimized variable of the function [x1,x2,x3,....].
        Lower,Upper : array
            The lower and upper limit of the variable.
        dcon : array, optional
            The derivatives of the constrains. The default is dcon.
        m : int, optional
            Number of constrains. The default is vel.
        c0 : array , optional
            intial values of the compliance. The default is c0.
        Xmin,Xmax : array, optional
            The maximuum and minimum values of X. The default is Xmin.
        tol : float, optional
             tolerance for minimizing of constrains. The default is 1e-5.
        tol1,tol2  : float, optional
            Tolerance of objective constrain. The default is 0.001 and 1.001.
        
        Returns
        -------
        p : array
            Function approximation for constrains equ. 4.2 in MMA paper.
        q : array
            Function approximation for constrains equ. 4.2 in MMA paper.
        rc: array
            residual of the function based on equ. 4.5 in MMA paper

        '''
        df_postive=np.maximum(dcon,0)  #maximum of the constrains based on equ. 4.3 in MMA paper
        df_negative=np.maximum(-dcon,0) #minimum of the constrains based on equ. 4.3 in MMA paper
        
        Xd_inv=np.array([1/Xmax-Xmin]).T
        UX=np.array([Upper-x]).T 
        LX=np.array([x-Lower]).T
        
        tol_vector=tol*np.ones((m,1)).T
        tol_matrix=np.dot(Xd_inv,tol_vector).T

        
        p0=tol2*df_postive+tol1*df_negative+tol_matrix
        
        UX_matrix=np.diag(((Upper-x)**2),0)
        LX_matrix=np.diag(((x-Lower)**2),0)
        
        p=(UX_matrix@p0.T).T  #Based on equ, 4.3 in MMA paper
        #print('p',p)
        
        q0=tol1*df_postive+tol2*df_negative+tol_matrix
        q=(LX_matrix@q0.T).T #Based on equ, 4.4 in MMA paper
        #print('q',q)
        
        pr=np.dot(p,(1/UX)) 
        qr=np.dot(q,(1/LX))
        rc= c0-(pr+qr).T ##Based on equ, 4.5 in MMA paper
        #print('rc',rc.T)
        return p,q,rc.T


    def prime_dual(L,U,alpha,beta,p0,q0,pc,qc,a,b,c,d,n=nel,m=vel,epsimin=1e-7):
        '''
        This function is a active set strategy and solves constrained based non-linear problem. It is a gradient based method using Newton Raphson to solve KKT condition.

        Parameters
        ----------
        L,U : array
            Lower and upper limits of the variables.
        alpha,beta : array
            Bound of the variables . Based on equ. 3.6 and 3.7 in MMA paper
        p0,q0 : array
            The approximation function of objectives based on equ. 3.3 and 3.4.
        pc,qc : array
            The approximation function of constraind based on equ. 3.3 and 3.4.
        a,b,c,d : array
            Non-negative lagragian multiplier.
        n : int, optional
            Number of variables. The default is nel.
        m : int, optional
            Number of constrains. The default is vel.
        epsimin : float, optional
            Tolerance value. The default is 1e-7.

        Returns
        -------
        x : array
            The value of variable obtained after primal dual method
        ii1: int
            Number of iterations taken by newton method to solve the problem
        l_ii1: int
            Total number of iterations taken by newton method to solve the problem
        res : array
            residual obtained from solving KKT conditions.

        '''
        def initial_condition(alpha=alpha,beta=beta,o=n,k=m,c=c):
            '''
            For solving constrained based problem. The initial values are required 

            Parameters
            ----------
            alpha,beta :array
               Bound of the variables . Based on equ. 3.6 and 3.7 in MMA paper
            o : int
                number of variables 
            k : int
                Number of constrains
            c : array, optional
                Non-negative Lagrangian multiplier

            Returns
            -------
            x,y,z : array
                The values of lagrangian function variables based on equ. 5.6 in MMA paper

            epsi,lamda,s : array
                weight parameters in lagrangian function.
            Zee,eta,neta,nu : array
                Non-negative lagrangian multiplier in equ. 5.6 in MMA paper.
            variable_vector,minimizer_vector : array
                initial value vector for weights(variables) and Lagrangian multilpier(constrains) initialized with 1.
            '''
            #the initial condition are obtained from MMA paper section 5.1
            x=np.array([0.5*(alpha+beta)]).T  #for variables
            y=np.ones((k,1)) # for constrains
            z=np.array([[1.0]])
            #print(o,k)
            #initial values for weight and Lagrangian multipler
            variable_vector=np.ones((o,1)) # based on number of variables
            minimizer_vector=np.ones((k,1)) #based on number of constrains
            #weight parameter in equ. 5.6 in MMA paper
            epsi=1
            lamda=minimizer_vector
            s=minimizer_vector
            #Non-negative Lagrangian multipler 
            Zee=z
            eta=np.maximum(1,1/(x-np.array([alpha]).T))
            neta=np.maximum(1,1/(((np.array([beta]).T)-x)))
            nu=np.maximum(lamda,c*0.5)
            return x,y,z,epsi,lamda,s,Zee,eta,neta,nu,variable_vector,minimizer_vector
    
        x,y,z,epsi,lamda,s,Zeta,eta,neta,nu,variable_vector,minimizer_vector=initial_condition() #intial condition are initialized
        #print('initial_condition',initial_condition())
        def optimal_condtition(x,y,z,lamda,s,Zeta,eta,neta,nu,epsi,o='o',U=U,L=L,p0=p0,q0=q0,pc=pc,qc=qc,b=b,variable_vector=variable_vector,minimizer_vector=minimizer_vector):
            '''
            this function calculates KKT conditions.

            Parameters
            ----------
            x,y,z : array
                 The values of lagrangian function variables based on equ. 5.6 in MMA paper
            epsi,lamda,s : array
                weight parameters in lagrangian function.
            Zee,eta,neta,nu : array
                Non-negative lagrangian multiplier in equ. 5.6 in MMA paper.
            L,U : array
            Lower and upper limits of the variables.
            p0,q0 : array
                The approximation function of objectives based on equ. 3.3 and 3.4.
            pc,qc : array
                The approximation function of constraind based on equ. 3.3 and 3.4.
            b: array
                    residual of approximate function equ. 3.2 in MMA paper.
            variable_vector,minimizer_vector : array
                initial value vector for weights(variables) and Lagrangian multilpier(constrains) initialized with 1.
            Returns
            -------
            residumax : array
                The maximum residual obtained after solving conditions 5.7 in MMA paper.
            residunorm : array
                Norm of the residuals obtained after solving conditions 5.7 in MMA paper.

            '''
            # based on equ. 5.5 in MMA paper
            UX=np.array([U]).T-x
            LX=x-np.array([L]).T
            #print(pc.T*lamda)
            plamda=p0+np.dot(pc.T,lamda)
                #print('plamda',plamda)
            qlamda=q0+np.dot(qc.T,lamda)
            g=np.dot(pc,(1/(UX)))+np.dot(qc,(1/(LX)))
            dphi_dx=(plamda/UX**2)-(qlamda/LX**2) #Based on equ. 5.4 in MMA paper
            #print('g',g)
            #print('dphi_dx',dphi_dx)
            #optimality conditions 5.7
            dl_dx=dphi_dx-eta+neta
            #print('dl_dx',dl_dx)
            dl_dy=c+d*y-lamda-nu
            dl_dz=a0-Zeta-np.dot(a.T,lamda)
            #print('dl_dy',dl_dy)
            #print('dl_dz',dl_dz)
            residu1 = np.concatenate((dl_dx, dl_dy, dl_dz),axis=0)
            #print('residu1',residu1)
            
            
            #residuals calculation
            #print('a',a)
            #print('z',z)
            #print(np.dot(a.T,z))

            #residual are calculated based on equ. 5.7 in MMA paper
            r1=g-np.dot(a,z)-y+s-b
            #print('r1',r1)
            r2=(eta*(x-np.array([alpha]).T))-epsi*variable_vector 
            #print('r2',r2)
            r3=(neta*(np.array([beta]).T-x))-epsi*variable_vector
            #print('r3',r3)
            r4=nu*y-epsi*minimizer_vector
            #print('r4',r4)
            r5=Zeta*z-epsi
            #print('r5',r5)
            r6=lamda*s-epsi*minimizer_vector
            #print('r6',r6)
            #print(r1,r2,r3)
            residu2 = np.concatenate((r1, r2, r3, r4, r5, r6), axis = 0)
            residu = np.concatenate((residu1, residu2), axis = 0)
            #print('residu',residu.T)
            residunorm = np.sqrt((np.dot(residu.T,residu)).item())
            #print(residunorm)
            residumax = np.max(np.abs(residu))
            #print(residumax)
            return residumax,residunorm
        
        def line_search(w,dw,x,y,z,lamda,eta,neta,nu,Zeta,s,dx,dy,dz,dlamda,deta,dneta,dnu,dZeta,ds,resi_norm,epsi,p0=p0,q0=q0,pc=pc,qc=qc,U=U,L=L,alpha=alpha,beta=beta,variable_vector=variable_vector,minimizer_vector=minimizer_vector):
            '''
            Newton line search algorthim built based on setion 5.4 in MMA paper. Used to solve Newton raphson but find the correct step size and calculating KKT condition 

            Parameters
            ----------
            w,dw : array
                The values of variables like weight parameters and non-negative lagrangian multipltier and it changes calculated from solving linear equations.
            x,y,z : array
                 The values of lagrangian function variables based on equ. 5.6 in MMA paper
            epsi,lamda,s : array
                weight parameters in lagrangian function.
            Zee,eta,neta,nu : array
                Non-negative lagrangian multiplier in equ. 5.6 in MMA paper.
            L,U : array
            Lower and upper limits of the variables.
            p0,q0 : array
                The approximation function of objectives based on equ. 3.3 and 3.4.
            pc,qc : array
                The approximation function of constraind based on equ. 3.3 and 3.4.
            b: array
                    residual of approximate function equ. 3.2 in MMA paper.
            variable_vector,minimizer_vector : array
                initial value vector for weights(variables) and Lagrangian multilpier(constrains) initialized with 1.
            dx,dy,dz,dlamda,dneta,dzeta,ds : array
                obtained from equ. 5.12 and 5.13 in MMA paper.
            resi_norm : float
                Previous iteration's residual.
            epsi : float 
                maximum Tolerance.

            Returns
            -------
            newx,newy,newlam,neweta,newnu,new_s,newZeta : TYPE
                Updated values after newton line search which satisfy KKT conditions.

            new_residualmax : int
                Max of the residuals obtained from optimal conditions function.
            new_residunorm : int
                Norm of residuals.
            it : int
                Number of iteration taken by line search to obtain optimal conditions.

            '''
            #calculating initial step step
            #dw=np.maximum(dw,1e-5)
            #print(dw)
            #step size is caculated based on section 5.4 in MMA paper
            size = -1.01*dw/w

            max_size = np.max(size) 
            stepa =-1.01*(dx/(x-np.array([alpha]).T))
            max_stepa = np.max(stepa)
            stepb = 1.01*(dx/((np.array([beta]).T)-x))
            max_stepb = np.max(stepb)
            max_step = max(max_stepa,max_stepb)
            final_step = max(max_step,max_size)
            step_size = max(final_step,1.0)
            a = 1.0/step_size
            #print('a',a)
            it=0
            max_iteration=200
            while it<max_iteration:
                newx=x+a*dx
                newy=y+a*dy
                newz=z+a*dz
                newlam=lamda+a*dlamda
                neweta = eta+a*deta
                newneta = neta+a*dneta
                newnu = nu+a*dnu
                newZeta = Zeta+a*dZeta
                new_s = s+a*ds
            
                #print('lam,eta,neta,nu,Zeta,s',x,y,z,lamda,s,Zeta,eta,neta,nu)
                #KKT condition are checked and residuals are calculated
                new_residualmax,new_residunorm=optimal_condtition(newx,newy,newz,newlam,new_s,newZeta,neweta,newneta,newnu,epsi)
                it+=1
                
                
                #print('line_residual',[new_residualmax,new_residunorm])
                if new_residunorm <(2*resi_norm): #exit condition (if norm of residual is less then 2* initial norm of residual)
                      return newx,newy,newz,newlam,neweta,newneta,newnu,newZeta,new_s,new_residualmax,new_residunorm,it 
                a=a*0.5 #change in step size if the solution doesnot converge.

            
        
        def Newton_method(U,L,x,y,z,alpha,beta,p0,q0,pc,qc,epsi,lamda,s,Zeta,eta,neta,nu,residumax,residunorm,variable_vector=variable_vector,minimizer_vector=minimizer_vector):
            '''
            Gradient based method where newton line search is employed to find the optimal search direction and step
            A system of linear equation are build based on equ. 5.11, 5.12, 5.13 in MMA paper.

            Parameters
            ----------
            L,U : array
            Lower and upper limits of the variables.
            x,y,z : array
                 The values of lagrangian function variables based on equ. 5.6 in MMA paper
            alpha,beta : array
                Bound of the variables . Based on equ. 3.6 and 3.7 in MMA paper
            p0,q0 : array
                The approximation function of objectives based on equ. 3.3 and 3.4.
            pc,qc : array
                The approximation function of constraind based on equ. 3.3 and 3.4.
            epsi,lamda,s : array
                weight parameters in lagrangian function.
            Zee,eta,neta,nu : array
                Non-negative lagrangian multiplier in equ. 5.6 in MMA paper.
            residumax : array
                The maximum residual obtained after solving KKT conditions 5.7 in MMA paper.
            residunorm : array
                Norm of the residuals obtained after solving KKT conditions 5.7 in MMA paper.
            
            b: array
                    residual of approximate function equ 3.2 .
            variable_vector,minimizer_vector : array
                initial value vector for weights(variables) and Lagrangian multilpier(constrains) initialized with 1.
            Returns
            -------
            x,y,z,lamda,eta,neta,nu,Zeta,s: array
                Updated values of variable, weights and lagrangian multiplier obtained after newton raphson
            iteration: int
                Number of iteration taken by Newton_raphson to achieve optimal conditons
            l_ii: 
              Total Number of iteration taken by newton line search to achieve optimal conditons
            
            '''
            def linear_system_assembly(Dx,Dy,Dlamda,delx,dely,delz,dellamda,G,a,z,Zeta,variable_vector=variable_vector,minimizer_vector=minimizer_vector):
                '''
                This function builts a system of linear equation based on matrix 5.14 

                Parameters
                ----------
                Dx,Dy,Dlamda,delx,dely,delz,dellamda : TYPE
                    Values required to build linear equation based on equ. 5.14.
                G,a,z,Zeta : array
                    Values required to build linear equation based on equ. 5.12,5.11
                variable_vector,minimizer_vector : array
                    initial value vector for weights(variables) and Lagrangian multilpier(constrains) initialized with 1.

                Returns
                -------
                dx,dy,dz,dlamda : array
                    The solution of system of linear equation.
                '''
                Dx_inv=1/Dx
                Dy_inv=1/Dy
                Dlamda_y=Dlamda+Dy_inv
                dellamda_y=dellamda+Dy_inv*dely
                #Building system of linear equation 5.20 in MMA paper
                if len(variable_vector)>len(minimizer_vector):
                    ###
                    A11=np.asarray(np.diag(Dlamda_y.flatten(),0) \
                        +(np.diag(Dx_inv.flatten(),0).dot(G.T).T).dot(G.T))
                    ###
                    #print(A11)
                    A12=a
                    A21=A12
                    A22=-Zeta/z
                    
                    A=np.concatenate((np.concatenate((A11,A12),axis=1),np.concatenate((A21,A22),axis=0).T),axis=0)
                    #print('A',A)
                    B1=dellamda_y-np.dot(G,(delx*Dx_inv))
                    B2=delz
                
                    B=np.concatenate((B1,B2),axis=0)
                    #print('B',B)
                    #solving system of eqautions
                    X=np.linalg.solve(A,B)
                    #obtaining results
                    dlamda=X[:len(minimizer_vector)]
                    dz=X[len(minimizer_vector):len(minimizer_vector)+1]
                    dx=-(Dx_inv*np.dot(G.T,dlamda))-Dx_inv*delx
                    dy=(Dy_inv*dlamda)-(Dy_inv*dely)
                    
                    return dx,dy,dz,dlamda
            #start of Newton Raphson method        
            residnorm=residunorm  #residual from previous iteration      
            iteration=0
            l_ii=0
            max_iteration=250 
            #exit conditions
            while iteration<max_iteration:
                iteration+=1
                UX=np.array([U]).T-x
                LX=x-np.array([L]).T
                plamda=p0+np.dot(pc.T,lamda)
                #print('plamda',plamda)
                qlamda=q0+np.dot(qc.T,lamda)
                ## based on equ. 5.4
                g=np.dot(pc,(1/(UX)))+np.dot(qc,(1/(LX)))
                dphi_dx=(plamda/UX**2)-(qlamda/LX**2)
                dphi_dxdx=(2*plamda/UX**3)+(2*qlamda/LX**3)
                
                #print('phi',g)
                #print('dphi_dxdx',dphi_dxdx)
                #dphi_dxdx=2*plamda.T*(1/((U-x)**3))+2*qlamda.T*(1/((x-L)**3))
                #print('dphi_dxdx',dphi_dxdx)
                ##based on equ. 5.12   
                Dx=dphi_dxdx+(eta/(x-np.array([alpha]).T))+(neta/(np.array([beta]).T-x))
                Dy=d+nu/y
                Dlamda=s/lamda
                #print('Dx',Dx)
                #print('Dy',Dy)
                #print('Dlamda',Dlamda)
                #print('dphi_dx',dphi_dx)
                ##Based on equ. 5.13 in MMA paper
                delx=dphi_dx-(1/(x-np.array([alpha]).T))*epsi*variable_vector+(1/(np.array([beta]).T-x))*epsi*variable_vector
                dely=c+d*y-lamda-(epsi*minimizer_vector)/y
                delz=a0-np.dot(lamda.T,a)-epsi/z
                dellamda=g-a*z-y-b+epsi*minimizer_vector/lamda
                #print('delx',delx)
                #print('dely',dely)
                #print('delz',delz)
                #print('dellamda',dellamda)
                Dx_inv=1/Dx
                Dy_inv=1/Dy
                Dlamda_y=Dlamda+Dy_inv
                dellamda_y=dellamda+Dy_inv*dely
   
                inv_UX=1/(U-x.T)**2
                diag_pij=inv_UX[0,:]
                inv_LX=1/(x.T-L)**2
                diag_qij=inv_LX[0,:]
                pij=(np.diag(diag_pij,0).dot(pc.T)).T
                qij=(np.diag(diag_qij,0).dot(qc.T)).T
                G=pij-qij #based on equ. 5.11
                #print('G',G)
                ## A system of linear equation is built 5.14 and solution is obtained which is search direction
                dx,dy,dz,dlamda=linear_system_assembly(Dx,Dy,Dlamda,delx,dely,delz,dellamda,G,a,z,Zeta)
                #print('dx',dx,'dy',dy,'dz',dz,'dlamda',dlamda)
                ##Non-negative lagrangian multipler are calculated based on equ. 5.13
                AlphaX=x-np.array([alpha]).T
                BetaX=np.array([beta]).T-x
                deta=-((eta*dx)/AlphaX)-eta+((epsi*variable_vector)/AlphaX)
                dneta=((neta*dx)/BetaX)-neta+((epsi*variable_vector)/BetaX)
                #print('dneta',dneta)
                dnu=-(nu*dy/y)-nu+(epsi*minimizer_vector/y)
                dZeta=-((Zeta/z)*dz)-Zeta+(epsi/z)
                ds=-((s*dlamda)/lamda)-s+epsi*minimizer_vector/lamda
                ##The equation is built based on section 5.3 in MMA paper
                w=np.concatenate((y,z,lamda,eta,neta,nu,s,Zeta),axis=0)
                dw=np.concatenate((dy,dz,dlamda,deta,dneta,dnu,ds,dZeta),axis=0)
                ##Inputs for newton line search 
                oldx=x
                oldy=y
                oldz=z
                oldlamda=lamda
                oldeta=eta
                oldneta=neta
                oldnu=nu
                oldZeta=Zeta
                olds=s
                #Newton line search algorthim is built and search direction and pervious iteration values are the input
                new_x,new_y,new_z,new_lamda,new_eta,new_neta,new_nu,new_Zeta,new_s,new_residualmax,new_residunorm,ii2=line_search(w,dw,x,y,z,lamda,eta,neta,nu,Zeta,s,dx,dy,dz,dlamda,deta,dneta,dnu,dZeta,ds,residnorm,epsi)
                l_ii+=ii2
                #updating variables , weight parameters, lagrangian multipler,residuals from results obtained from newton linea search
                x=new_x
                y=new_y
                z=new_z
                lamda=new_lamda
                eta=new_eta
                neta=new_neta
                nu=new_nu
                Zeta=new_Zeta
                s=new_s
                residumax=new_residualmax
                residnorm=new_residunorm
                #calculating the change
                xx1=oldx-x
                xx2=oldy-y
                xx3=oldz-z
                xx4=oldlamda-lamda
                xx5=oldeta-eta
                xx6=oldneta-neta
                xx7=oldnu-nu
                xx8=oldZeta-Zeta
                xx9=olds-s

                exit_cond=np.concatenate((xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8,xx9),axis=0)
                if iteration>=max_iteration:
                    print('line_search',exit_cond.T)
                    print('w',w.T)
                    print('dw',dw.T)
                #print('line residumax:',[residumax,residnorm])
                ##Newton raphson would exit when the maximum residual is less than given tolerance 
                if residumax<0.9*epsi :#or np.linalg.norm(abs(exit_cond))<1e-4:

                    #print(x,y,z,lamda,eta,neta,nu,Zeta,s)
                    
                    return x,y,z,lamda,eta,neta,nu,Zeta,s,iteration,l_ii
        
        ii1=0
        l_ii1=0
        #start of prime dual method
        while epsi> epsimin: #exit condition
            #print('x',x)
            #KKT conditions are calculated
            residual_max,residunorm=optimal_condtition(x,y,z,lamda,s,Zeta,eta,neta,nu,epsi)
            #print('outerresidual:',[residual_max,residunorm])
            oldx=x
            oldy=y
            oldz=z
            oldlamda=lamda
            oldeta=eta
            oldneta=neta
            oldnu=nu
            oldZeta=Zeta
            olds=s
            #Newton-Raphson is performed
            x,y,z,lamda,eta,neta,nu,Zeta,s,ii1,l_ii=Newton_method(Upper,Lower,x,y,z,alpha,beta,p0,q0,pc,qc,epsi,lamda,s,Zeta,eta,neta,nu,residual_max,residunorm)
            ii1+=ii1
            l_ii1+=l_ii
            #change is calculated
            xx1=oldx-x
            xx2=oldy-y
            xx3=oldz-z
            xx4=oldlamda-lamda
            xx5=oldeta-eta
            xx6=oldneta-neta
            xx7=oldnu-nu
            xx8=oldZeta-Zeta
            xx9=olds-s
            epsi=epsi*0.1 
            exit_cond=np.concatenate((xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8,xx9),axis=0)
            res=np.linalg.norm(abs(exit_cond))
            #if np.linalg.norm(abs(exit_cond))<1e-4:
            #    print(x)
        return x,ii1,l_ii1,res
        
    Lower,Upper,x=Asymptoes()
    c = 10000*np.ones((vel,1))
    #print(c)
    d = np.ones((vel,1))
    a0 = 1
    a = np.zeros((vel,1))
    #calculating alpha and beta
    alpha=np.maximum(Xmin,np.maximum((Lower+m_tol*(x-Lower)),(x-m*(Xmax-Xmin))))
    beta=np.minimum(Xmax,np.minimum((Upper-m_tol*(Upper-x)),(x+m*(Xmax-Xmin))))
    #calculating dervivative of objective and constrains 
    p0,q0=objective_constrains(x,Lower,Upper)
    pc,qc,rc=Minimizer_constrains(x,Lower,Upper)
    b=-rc #residual of approximate function
    #print(Lower,Upper,alpha,beta,p0,q0,pc,qc,a,b,c,d)
    ##constrained based problem is solved using a gradient based active set strategy method 
    x,exit_i,l_iit,res=prime_dual(Lower,Upper,alpha,beta,p0,q0,pc,qc,a,b,c,d)
    #print(x)
    
    TGREEN = '\033[32;1m'
    TBLUE = '\033[34;1m'
    ENDC = '\033[m'
    if True:
        width=120
        print('*'*width)
        fmt='{:^'+str(width)+'}'
        print(TBLUE+fmt.format('Method of Moving Asymptoes')+ENDC)
        print('`'*width)
        fmt='{:^'+str(width)+'}'
        print(TGREEN+fmt.format('Prime-Dual Newton Method')+ENDC)
        print('-'*width)
        print(TBLUE+'             Lower: %f3          Upper:%f3           Alpha:%f3           Beta:%f3 ' %(np.mean(Lower),np.mean(Upper),np.mean(alpha),np.mean(beta))+ENDC) 
        print('-'*width)
        print(TGREEN+'     Newton_iter: %d      Line_search_ite:%d      residual:%f4      vol_constrain:%f6     vol:%f4' %(exit_i,l_iit,res,0,np.mean(x))+ENDC)  
        print('*'*width)
        print('#'*width)
    return x,Lower,Upper
