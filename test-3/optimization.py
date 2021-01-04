import numpy as np
import pytest

def Knearestneighbours(rmin,nelx,nely,nelz): ##modify
    nel=nelx*nely*nelz
    H=np.zeros((nel,nel))
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
                            H[r,c]=np.maximum(0,(rmin-distance))
    H=np.array(H)
    DH=np.sum(H,1)
    return H,DH  


def optimality_criteria(dfun,dCon,initial_X,constrain,H,DH,beta,oc_disp=True,g=1,tol=0.01,move=0.1,neta=0.5,initial_value=0,final_value=1e09):
    X_new=np.zeros(len(initial_X))
    #X_filter=np.zeros(len(initial_X))
    X=np.array(initial_X)
    i=0
    max_iteration=150
    while i<=max_iteration:
        convergence_criteria=(final_value-initial_value)/(final_value+initial_value)

        #lagarian multiplier
        lamdba=0.5*(initial_value+final_value)
        #Bey
        Be=((-dfun/(dCon*lamdba))**neta)**g
        X_new[:]= np.maximum(0.0,np.maximum(X-move,np.minimum(1.0,np.minimum(X+move,X*Be))))
        #X_new=np.matmul(H,X_new/DH)
        #X_new=1-np.exp(-beta*X_filter)+X_filter*np.exp(-beta)
        update_condition=np.sum(X_new)
        #  bi-section algorthim updating lagaranian multiplier
        if update_condition>=constrain:
            initial_value=lamdba    
        else:
            final_value=lamdba
        i=i+1
        if convergence_criteria<=tol :#or (np.linalg.norm(abs(X_new-X_old)))>=1e-4:
            width=120
            if oc_disp:
                print('-'*width)
                fmt='{:^'+str(width)+'}'
                print(fmt.format('Optimiality Criterian'))
                print('.'*width)
                print('     exit_iter: %d    gray_filter(g):%f    lamdba:%f12      tot_vol:%f4     max_volume:%f4' %(i,g,lamdba,update_condition,constrain))  
                print('.'*width)
            return X_new
        if i==(max_iteration-1):
            width=120
            if oc_disp:
                print('maximum iteration has reached, didnot converger')
                print('.'*width)
                fmt='{:^'+str(width)+'}'
                print(fmt.format('Optimiality Criterian'))
                print('.'*width)
                print('     exit_iter: %d        lamdba:%f12        tot_vol:%f4       max_volume:%f4' %(i,lamdba,update_condition,constrain))  
                print('*'*width)
            return X_new

#def Moving_asymptoes

def Heavyside_filter(element_density,beta,ll=0,lh=1,tol=0.01,fil_disp=True):
    be_element=np.zeros(len(element_density))
    i=0
    max_iteration=100
    while i<max_iteration:
        cv=(ll+lh)/2
        be_element=(np.tanh((beta*cv))+np.tanh(beta*(element_density-cv))/(np.tanh((beta*cv))+np.tanh(beta*(1-cv))))
        filtered_density=np.maximum(1e-3,be_element)
        if sum(be_element)>sum(element_density):
            ll=cv
        else:
            lh=cv
        i+=1
        #if abs((lh-ll))<1e-5:
        if abs((lh-ll))<1e-5:
            if fil_disp:
                width=120
                fmt='{:^'+str(width)+'}'
                print(fmt.format('Heavyside Filter'))
                print('*'*width)
                print('         exit_iter: %d          n:%f12           vol:%f4           max_volume:%f4' %(i,cv,sum(be_element),sum(element_density)))  
                print('.'*width)
                return filtered_density

#element_density=np.array([1,1,0,1,.5,1,0.2,1])
#element_density=np.random.random((10,1))
#beta=1
#print(Heavyside_filter(element_density,beta,ll=0,lh=1,tol=0.01,fil_disp=True))
#print(element_density)


def Moving_asymptoes(dfun0,dcon,f0,c0,x0,x1,x2,L,U,loop,nel,vel,Xmin,Xmax,ma_disp,m=0.2,m_tol=0.1):
    def Asymptoes(loop=loop,x0=x0,x1=x1,x2=x2,n=nel,Xmin=Xmin,Xmax=Xmax,L=L,U=U,move_tol_low=0.01,move_tol_high=10):#
        Lower=np.ones(n)
        Upper=np.ones(n)
        if loop<=2:
            Lower=x0-0.5*(Xmax-Xmin)
            Upper=x0+0.5*(Xmax-Xmin)
            #print(Lower,Upper,x)
        else:
            '''
            gamma=(x0-x1)*(x1-x2)
            gamma[gamma<0]=0.7
            gamma[gamma>0]=1.2
            #gamma[gamma==0]=0

            Lower=x0-0.5*gamma*(x1-L)
            Upper=x0+0.5*gamma*(U-x2)
        
            Lower=np.maximum(Lower,x0-move_tol_high)
            Lower=np.minimum(Lower,x0-move_tol_low)
            Upper=np.maximum(Upper,x0+move_tol_low)
            Upper=np.minimum(Upper,x0+move_tol_high)
            '''
            xval=x0
            xold1=x1
            xold2=x2
            zzz = (xval-xold1)*(xold1-xold2)
            gamma = np.ones(nel)
            gamma[np.where(zzz>0)] = 1.2
            gamma[np.where(zzz<0)] = 0.7
            #gamma[np.where(zzz==0)]= 0 
            Lower = xval-gamma*(xold1-L)
            #print(Lower)
            Upper = xval+gamma*(U-xold1)
            #print(Upper)
            lowmin = xval-10*(Xmax-Xmin)
            lowmax = xval-0.01*(Xmax-Xmin)
            uppmin = xval+0.01*(Xmax-Xmin)
            uppmax = xval+10*(Xmax-Xmin)
            Lower = np.maximum(Lower,lowmin)
            Lower = np.minimum(Lower,lowmax)
            Upper = np.minimum(Upper,uppmax)
            Upper = np.maximum(Upper,uppmin)
        ##output as : [ 1.5  0.5 -0.5], [6.5 5.5 4.5], [4 3 2]
        return Lower,Upper,x0
    
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

    def Minimizer_constrains(x,Lower,Upper,dcon=dcon,m=vel,c0=c0,Xmin=Xmin,Xmax=Xmax,tol=1e-5,tol1=0.001,tol2=1.001):#
        df_postive=np.maximum(dcon,0)
        df_negative=np.maximum(-dcon,0)
        
        Xd_inv=np.array([1/Xmax-Xmin]).T
        UX=np.array([Upper-x]).T
        LX=np.array([x-Lower]).T
        
        tol_vector=tol*np.ones((m,1)).T
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


    def prime_dual(L,U,alpha,beta,p0,q0,pc,qc,a,b,c,d,n=nel,m=vel,epsimin=1e-7):
        
        def initial_condition(alpha=alpha,beta=beta,o=n,k=m,c=c):
            x=np.array([0.5*(alpha+beta)]).T
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
        def optimal_condtition(x,y,z,lamda,s,Zeta,eta,neta,nu,epsi,o='o',U=U,L=L,p0=p0,q0=q0,pc=pc,qc=qc,b=b,variable_vector=variable_vector,minimizer_vector=minimizer_vector):
            #5.5
            UX=np.array([U]).T-x
            LX=x-np.array([L]).T
            #print(pc.T*lamda)
            plamda=p0+np.dot(pc.T,lamda)
                #print('plamda',plamda)
            qlamda=q0+np.dot(qc.T,lamda)
                #5.4
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
            #calculating initial step step
            #dw=np.maximum(dw,1e-5)
            #print(dw)
            '''
            step_size=np.max(-1.01*(dw/w))
            #print(step_size)
            DL=np.max(-1.01*(dx/(x-np.array([alpha]).T)))
            DU=np.max(1.01*(dx/(np.array([beta]).T)-x))
            step_range=max(DL,DU)
            a=1/max(1.0,max(step_size,step_range))
            '''
            stepxx = -1.01*dw/w
            #out = -1.01*np.ones( (8) )  #preinit
            #stepxx=-1.01*np.divide(dw, w, out=np.zeros_like(dw), where=w!=0)
            stmxx = np.max(stepxx) 
            stepalfa =-1.01*(dx/(x-np.array([alpha]).T))
            stmalfa = np.max(stepalfa)
            stepbeta = 1.01*(dx/((np.array([beta]).T)-x))
            stmbeta = np.max(stepbeta)
            stmalbe = max(stmalfa,stmbeta)
            stmalbexx = max(stmalbe,stmxx)
            stminv = max(stmalbexx,1.0)
            a = 1.0/stminv
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
                #np.set_printoptions(precision=4)
                #print('x,y,z',newx,newy,newz)
                #print('\n lamda',newlam)
                #print('\n neweta',neweta)
                #print('\n newneta',newneta)
                #print('\n newnu',newnu)
                #print('\n newZeta',newZeta)
                #print('\n new_s',new_s)
            
                #print('lam,eta,neta,nu,Zeta,s',x,y,z,lamda,s,Zeta,eta,neta,nu)
                new_residualmax,new_residunorm=optimal_condtition(newx,newy,newz,newlam,new_s,newZeta,neweta,newneta,newnu,epsi)
                it+=1
                
                
                #print('line_residual',[new_residualmax,new_residunorm])
                if new_residunorm <(2*resi_norm):
                      return newx,newy,newz,newlam,neweta,newneta,newnu,newZeta,new_s,new_residualmax,new_residunorm,it 
                a=a*0.5
                #print('a1',a)
            
        
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
            l_ii=0
            max_iteration=250
            while iteration<max_iteration:
                iteration+=1
                UX=np.array([U]).T-x
                LX=x-np.array([L]).T
                #print('x',x)
                #print('\n U',np.array([U]).T)
                #print('\n UX',UX)
                #print('\n UX',UX**3)
                #print('\n UX_INV',1/UX)
                #print('\n UX_INV2',(1/(UX**2)))
                #print('\n UX_INV3',(1/(UX**3)))
                plamda=p0+np.dot(pc.T,lamda)
                #print('plamda',plamda)
                qlamda=q0+np.dot(qc.T,lamda)
                #5.4
                g=np.dot(pc,(1/(UX)))+np.dot(qc,(1/(LX)))
                dphi_dx=(plamda/UX**2)-(qlamda/LX**2)
                dphi_dxdx=(2*plamda/UX**3)+(2*qlamda/LX**3)
                
                #print('phi',g)
                #print('dphi_dxdx',dphi_dxdx)
                #dphi_dxdx=2*plamda.T*(1/((U-x)**3))+2*qlamda.T*(1/((x-L)**3))
                #print('dphi_dxdx',dphi_dxdx)
                   
                Dx=dphi_dxdx+(eta/(x-np.array([alpha]).T))+(neta/(np.array([beta]).T-x))
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
                #print('G',G)
                
                dx,dy,dz,dlamda=linear_system_assembly(Dx,Dy,Dlamda,delx,dely,delz,dellamda,G,a,z,Zeta)
                #print('dx',dx,'dy',dy,'dz',dz,'dlamda',dlamda)
                AlphaX=x-np.array([alpha]).T
                BetaX=np.array([beta]).T-x
                deta=-((eta*dx)/AlphaX)-eta+((epsi*variable_vector)/AlphaX)
                dneta=((neta*dx)/BetaX)-neta+((epsi*variable_vector)/BetaX)
                #print('dneta',dneta)
                dnu=-(nu*dy/y)-nu+(epsi*minimizer_vector/y)
                dZeta=-((Zeta/z)*dz)-Zeta+(epsi/z)
                ds=-((s*dlamda)/lamda)-s+epsi*minimizer_vector/lamda
                
                w=np.concatenate((y,z,lamda,eta,neta,nu,s,Zeta),axis=0)
                dw=np.concatenate((dy,dz,dlamda,deta,dneta,dnu,ds,dZeta),axis=0)
                #lline search 
                #print('w',np.round(w.T))
                #print('dw',np.round(dw.T))
                oldx=x
                oldy=y
                oldz=z
                oldlamda=lamda
                oldeta=eta
                oldneta=neta
                oldnu=nu
                oldZeta=Zeta
                olds=s
                new_x,new_y,new_z,new_lamda,new_eta,new_neta,new_nu,new_Zeta,new_s,new_residualmax,new_residunorm,ii2=line_search(w,dw,x,y,z,lamda,eta,neta,nu,Zeta,s,dx,dy,dz,dlamda,deta,dneta,dnu,dZeta,ds,residnorm,epsi)
                l_ii+=ii2
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
                if residumax<0.9*epsi :#or np.linalg.norm(abs(exit_cond))<1e-4:
                    #print(x,y,z,lamda,eta,neta,nu,Zeta,s)
                    
                    return x,y,z,lamda,eta,neta,nu,Zeta,s,iteration,l_ii
        
                
               

        ii1=0
        l_ii1=0
        while epsi> epsimin:
            #print('x',x)
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
            x,y,z,lamda,eta,neta,nu,Zeta,s,ii1,l_ii=Newton_method(Upper,Lower,x,y,z,alpha,beta,p0,q0,pc,qc,epsi,lamda,s,Zeta,eta,neta,nu,residual_max,residunorm)
            ii1+=ii1
            l_ii1+=l_ii
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
    #x=Newton_Method(dfun0,f0,x,Lower,Upper,alpha,beta)
    #calculating dervivative of objective and constrains 
    p0,q0=objective_constrains(x,Lower,Upper)
    pc,qc,rc=Minimizer_constrains(x,Lower,Upper)
    b=-rc
    #print(L,U)
    #print(p0,q0,r0)
    #print(pc,qc,b)
    #print(a,b,c,d)
    #print(Lower,Upper,alpha,beta,p0,q0,pc,qc,a,b,c,d)
    #x,y,z,lamda,eta,neta,nu,Zeta,s=subsolv(2,3,L,U,alpha,beta,p0,q0,pc,qc,a0,a,b,c,d)
    x,exit_i,l_iit,res=prime_dual(Lower,Upper,alpha,beta,p0,q0,pc,qc,a,b,c,d)
    #x,y,z,lamda,eta,neta,nu,Zeta,s,residual=prime_dual(L,U,alpha,beta,p0,q0,pc,qc,a,b,c,d)
    #print(x)
    if True:
        width=120
        print('*'*width)
        fmt='{:^'+str(width)+'}'
        print(fmt.format('Method of Moving Asymptoes'))
        print('`'*width)
        fmt='{:^'+str(width)+'}'
        print(fmt.format('Prime-Dual Newton Method'))
        print('-'*width)
        print('             Lower: %f3          Upper:%f3           Alpha:%f3           Beta:%f3 ' %(np.mean(Lower),np.mean(Upper),np.mean(alpha),np.mean(beta))) 
        print('-'*width)
        print('     Newton_iter: %d      Line_search_ite:%d      residual:%f4      vol_constrain:%f6     vol:%f4' %(exit_i,l_iit,res,c0,np.mean(x)))  
        print('*'*width)
        print('#'*width)
    return x,Lower,Upper
