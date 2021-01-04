import numpy as np
def Moving_asymptoes(dfun0,dcon,f0,c0,x0,L,U,loop,nel,ma_disp,a0,a,d,c,Xmin,Xmax,move_tol_low=0.01,move_tol_high=10,m=0.5,m_tol=0.1):

    def Asymptoes(loop=loop,x0=x0,n=nel,Xmin=Xmin,Xmax=Xmax,L=L,U=U):#
        Lower=np.ones(n)
        Upper=np.ones(n)
        if loop<3:
            x=x0[:,loop]
            Lower=x-0.5*(Xmax-Xmin)
            Upper=x+0.5*(Xmax-Xmin)
            #print(Lower,Upper,x)
        else:
            xk_2=x0[:,loop-2]
            xk_1=x0[:,loop-1]
            x=x0[:,loop]
            gamma=(x-xk_1)*(xk_1-xk_2)
            gamma[gamma<0]=0.7
            gamma[gamma>0]=1.2
            gamma[gamma==0]=0

            Lower=x-0.5*gamma*(xk_1-L)
            Upper=x+0.5*gamma*(U-xk_2)
        
            Lower=np.maximum(Lower,x-move_tol_high)
            Lower=np.minimum(Lower,x-move_tol_low)
            Upper=np.maximum(Upper,x+move_tol_low)
            Upper=np.minimum(Upper,x+move_tol_high)
        ##output as : [ 1.5  0.5 -0.5], [6.5 5.5 4.5], [4 3 2]
        return Lower,Upper,x
    
    def objective_constrains(x,Lower,Upper,dfun0=dfun0,f0=f0,Xmin=Xmin,Xmax=Xmax,tol=1e-5,tol1=0.001,tol2=1.001):#
        df_postive=np.maximum(dfun0,0)
        df_negative=np.maximum(-dfun0,0)
        #print(df_postive)
        Xd_inv=1/Xmax-Xmin
        UX=Upper-x
        
        LX=x-Lower
        
        p0=tol2*df_postive+tol1*df_negative+(tol*Xd_inv)
        #print('p0',p0)
        p=np.array([(UX**2)*p0]).T
        q0=tol1*df_postive+tol2*df_negative+(tol*Xd_inv)
        q=np.array([(LX**2)*q0]).T

        #print('p',p)
        #print('q',q)    
        return p,q

    def Minimizer_constrains(x,Lower,Upper,dcon=dcon,c0=c0,Xmin=Xmin,Xmax=Xmax,tol=1e-5,tol1=0.001,tol2=1.001):#
        df_postive=np.maximum(dcon,0)
        df_negative=np.maximum(-dcon,0)
        
        Xd_inv=np.array([1/Xmax-Xmin]).T
        UX=np.array([Upper-x]).T
        LX=np.array([x-Lower]).T
        
        tol_vector=tol*np.array([np.ones((len(dcon)))])
        tol_matrix=np.dot(Xd_inv,tol_vector).T

        
        p0=tol2*df_postive+tol1*df_negative+tol_matrix
        
        UX_matrix=np.diag(((Upper-x)**2),0)
        LX_matrix=np.diag(((x-Lower)**2),0)
        
        p=(UX_matrix@p0.T).T
        #print('p',p)
        
        q0=tol1*df_postive+tol2*df_negative+tol_matrix
        q=(LX_matrix@q0.T).T
        #print('q',q)
        
        pr=np.dot(p,(1/UX))
        qr=np.dot(q,(1/LX))
        rc= c0-(pr+qr).T #S remove .T
        #print('rc',rc.T)
        return p,q,rc.T


    def prime_dual(L,U,alpha,beta,p0,q0,pc,qc,a,b,c,d,epsimin=1e-7):
        
        def initial_condition(alpha=alpha,beta=beta,rc=b,c=c):
            x=np.array([0.5*(alpha+beta)]).T
            o=len(alpha)
            k=len(rc)
            y=np.ones((k,1))
            z=np.array([[1.0]])
            #print(o,k)
            variable_vector=np.ones((o,1))
            minimizer_vector=np.ones((k,1))
            epsi=1
            lamda=minimizer_vector
            s=minimizer_vector
            Zee=z
            eta=np.maximum(1,1/(x-np.array([alpha]).T))
            neta=np.maximum(1,1/(((np.array([beta]).T)-x)))
            nu=np.maximum(lamda,c*0.5)
            return x,y,z,epsi,lamda,s,Zee,eta,neta,nu,variable_vector,minimizer_vector
    
        x,y,z,epsi,lamda,s,Zeta,eta,neta,nu,variable_vector,minimizer_vector=initial_condition()
        #print('initial_condition',initial_condition())
        def optimal_condtition(x,y,z,lamda,s,Zeta,eta,neta,nu,U=U,L=L,p0=p0,q0=q0,pc=pc,qc=qc,epsi=epsi,variable_vector=variable_vector,minimizer_vector=minimizer_vector):
            #5.5
            UX=np.array([U]).T-x
            LX=x-np.array([L]).T
            #print('UX',UX)
            #print('LX',LX)
            #print('p0',p0)
            #print('pc',pc)
            plamda=p0+np.dot(pc.T,lamda)
            #print('plamda',plamda)
            qlamda=q0+np.dot(qc.T,lamda)
            #5.
            #print(pc)
            #print(1/UX)
            #print(np.dot(pc,(1/(UX)).T))
            g=np.dot(pc,(1/(UX)))+np.dot(qc,(1/(LX)))
            dphi_dx=(plamda/UX**2)-(qlamda/LX**2)
            #print('g',g)
            #print('dphi_dx',dphi_dx)
            #optimality conditions 5.7
            #print('eta',eta)
            #print('neta',neta)
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
            r1=g-(a*z)-y+s-rc
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
            #print(residu)
            residunorm = np.sqrt((np.dot(residu.T,residu)).item())
            #print(residunorm)
            residumax = np.max(np.abs(residu))
            return residumax,residunorm
        
        def line_search(w,dw,x,y,z,lamda,eta,neta,nu,Zeta,s,dx,dy,dz,dlamda,deta,dneta,dnu,dZeta,ds,resi_norm,p0=p0,q0=q0,pc=pc,qc=qc,U=U,L=L,alpha=alpha,beta=beta,epsi=epsi,variable_vector=variable_vector,minimizer_vector=minimizer_vector):
            #calculating initial step step
            
            step_size=np.max(-1.01*(w/dw))
            DL=np.max(-1.01*(dx/(x-np.array([alpha]).T)))
            DU=np.max(1.01*(dx/(np.array([beta]).T)-x))
            step_range=max(DL,DU)
            a=1/max(1.0,max(step_size,step_range))
            it=0
            max_iteration=50

            while it<max_iteration:
                newx=x+a*dx
                newy=y+a*dy
                newz=z+a*dz
                newlam=lamda+a*dlamda
                neweta = eta+a*deta
                newneta = dneta+a*dneta
                newnu = nu+a*dnu
                newZeta = Zeta+a*dZeta
                new_s = s+a*ds
                new_residualmax,new_residunorm=optimal_condtition(newx,newy,newz,newlam,new_s,newZeta,neweta,newneta,newnu)
                
                if new_residunorm <(2*resi_norm):
                      return newx,newy,newz,newlam,neweta,newneta,newnu,newZeta,new_s,new_residualmax,new_residunorm 
                a=a*0.5
            
        
        def Newton_method(U,L,x,y,z,alpha,beta,p0,q0,pc,qc,epsi,lamda,s,Zeta,eta,neta,nu,residumax,residunorm,variable_vector=variable_vector,minimizer_vector=minimizer_vector):
            
            def linear_system_assembly(Dx,Dy,Dlamda,delx,dely,delz,dellamda,G,a,z,Zeta,variable_vector=variable_vector,minimizer_vector=minimizer_vector):
                Dx_inv=1/Dx
                Dy_inv=1/Dy
                Dlamda_y=Dlamda+Dy_inv
                dellamda_y=dellamda+Dy_inv*dely
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
                    X=np.linalg.solve(A,B)
                    dlamda=X[:len(minimizer_vector)]
                    dz=X[len(minimizer_vector):len(minimizer_vector)+1]
                    dx=-(Dx_inv*np.dot(G.T,dlamda))-Dx_inv*delx
                    dy=(Dy_inv*dlamda)-(Dy_inv*dely)
                    
                    return dx,dy,dz,dlamda
            residnorm=residunorm        
            iteration=0
            max_iteration=150
            while iteration<max_iteration:
                iteration+=1
                UX=np.array([U]).T-x
                LX=x-np.array([L]).T

                plamda=p0+np.dot(pc.T,lamda)
                #print('plamda',plamda)
                qlamda=q0+np.dot(qc.T,lamda)
                #5.4
                g=np.dot(pc,(1/(UX)))+np.dot(qc,(1/(LX)))
                dphi_dx=(plamda/UX**2)-(qlamda/LX**2)
                dphi_dxdx=(plamda/UX**3)+(qlamda/LX**3)
                #print('phi',g)
                #print('dphi_dx',dphi_dx)
                #dphi_dxdx=2*plamda*(1/((U-x)**3))+2*qlamda*(1/((x-L)**3))
                #print((1/(UX**2))
                   
                Dx=2*dphi_dxdx+(eta/(x-np.array([alpha]).T))+(neta/(np.array([beta]).T-x))
                Dy=d+nu/y
                Dlamda=s/lamda
                #print('Dx',Dx)
                #print('Dy',Dy)
                #print('Dlamda',Dlamda)
                #print('dphi_dx',dphi_dx)
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
                G=pij-qij
                #print(G)
                
                dx,dy,dz,dlamda=linear_system_assembly(Dx,Dy,Dlamda,delx,dely,delz,dellamda,G,a,z,Zeta)
                #print('dx',dx,'dy',dy,'dz',dz,'dlamda',dlamda)
                AlphaX=x-np.array([alpha]).T
                BetaX=np.array([beta]).T-x
                deta=-((eta*dx)/AlphaX)-eta+((epsi*variable_vector)/AlphaX)
                dneta=((neta*dx)/BetaX)-neta+((epsi*variable_vector)/BetaX)
                dnu=-(nu*dy/y)-nu+(epsi*minimizer_vector/y)
                dZeta=-((Zeta/z)*dz)-Zeta+(epsi/z)
                ds=-((s*dlamda)/lamda)-s+epsi*minimizer_vector/lamda
                
                w=np.concatenate((y,z,lamda,eta,neta,nu,s,Zeta),axis=0)
                dw=np.concatenate((dy,dz,dlamda,deta,dneta,dnu,ds,dZeta),axis=0)
                #lline search 
                #print('w',w)
                #print('dw',dw)
                new_x,new_y,new_z,new_lamda,new_eta,new_neta,new_nu,new_Zeta,new_s,new_residualmax,new_residunorm=line_search(w,dw,x,y,z,lamda,eta,neta,nu,Zeta,s,dx,dy,dz,dlamda,deta,dneta,dnu,dZeta,ds,residnorm)
                print(new_y)
                residumax=new_residualmax
                residnorm=new_residunorm
                if residumax<0.9*epsi :
                    #print(x,y,z,lamda,eta,neta,nu,Zeta,s,residual)
                    return x,y,z,lamda,eta,neta,nu,Zeta,s,residual
                
               

        
        while epsi> epsimin:
            #print('p0',p0)
            residual_max,residunorm=optimal_condtition(x,y,z,lamda,s,Zeta,eta,neta,nu)
            print(residual_max,residunorm)
            x,y,z,lamda,eta,neta,nu,Zeta,s,residual=Newton_method(U,L,x,y,z,alpha,beta,p0,q0,pc,qc,epsi,lamda,s,Zeta,eta,neta,nu,residual_max,residunorm)
            epsi=epsi*0.1
        return x,y,z,lamda,eta,neta,nu,Zeta,s,residual
        
    L,U,x=Asymptoes()
    #calculating alpha and beta
    alpha=np.maximum(Xmin,np.maximum((L+m_tol*(x-L)),(x-m*(Xmax-Xmin))))
    beta=np.minimum(Xmax,np.minimum((U-m_tol*(U-x)),(x+m*(Xmax-Xmin))))
    #x=Newton_Method(dfun0,f0,x,Lower,Upper,alpha,beta)
    #calculating dervivative of objective and constrains 
    p0,q0=objective_constrains(x,L,U)
    pc,qc,rc=Minimizer_constrains(x,L,U)
    b=-rc
    #print(L,U)
    #print(p0,q0,r0)
    #print(pc,qc,b)
    #print(a,b,c,d)
    #print(L,U,alpha,beta,p0,q0,pc,qc,a,b,c,d)
    #x,y,z,lamda,eta,neta,nu,Zeta,s=subsolv(2,3,L,U,alpha,beta,p0,q0,pc,qc,a0,a,b,c,d)
    x,y,z,lamda,eta,neta,nu,Zeta,s,residual=prime_dual(L,U,alpha,beta,p0,q0,pc,qc,a,b,c,d)
    #x,y,z,lamda,eta,neta,nu,Zeta,s,residual=prime_dual(L,U,alpha,beta,p0,q0,pc,qc,a,b,c,d)

    return x,L,U


def quad_function(x):
        return x[0]**2+x[1]**2+x[2]**2

def dfunction(x):
    df_dx=np.array([2*x[0],2*x[1],2*x[2]])
    return df_dx
def constrain(x):
    c1=(x[0]-5)**2+(x[1]-2)**2+(x[2]-1)**2-9
    c2=(x[0]-3)**2+(x[1]-4)**2+(x[2]-3)**2-9
    C=np.array([c1,c2])
    return C
def dconstrains(x):
    c11=2*(x[0]-5)
    c12=2*(x[1]-2)
    c13=2*(x[2]-1)
    c21=2*(x[0]-3)
    c22=2*(x[1]-4)
    c23=2*(x[2]-3)
    dc=np.array([[c11,c12,c13],[c21,c22,c23]])
    return dc
x=np.array([4,3,2])
f0=quad_function(x)
df0=dfunction(x)
c0=constrain(x)
dc0=dconstrains(x)
#print(x,df0,dc0,f0,c0)


m = 2
n = 3
epsimin = 0.0000001
eeen = np.ones(n).T
eeem = np.ones(m).T
#print(eeem)
zeron = np.zeros(n)
zerom = np.zeros(m)
xval = np.array([[4,3,2]]).T
#print(xval)
xold1 = xval
xold2 = xval
xmin = zeron
xmax = 5*eeen
x_inv=np.array([1/(xmax-xmin)]).T

#print(x_inv)
L = xmin
U = xmax
move = 1.0
c = 1000*np.ones((m,1))
print(c)
d = np.ones((m,1))
a0 = 1
a = np.zeros((m,1))
#xmamiinv=np.array(1/(x-L))
#print(xmamiinv)
#print(np.dot(eeem.T,xmamiinv))
#p0 = p0*ux2
#q0 = q0*xl2
loop=0
nel=3

x,L,U=Moving_asymptoes(df0,dc0,f0,c0,xval,L,U,loop,nel,True,a0,a,d,c,xmin,xmax)
#print(x,L,U)