from numpy import linspace,array


def exact_displacements(P,l,b,h,E,v,nx,ny,nz):
    I=b*(h**3)/12
    disp=[]
    n=linspace(0,l,(nx))
    m=linspace(0,h,(ny))
    q=linspace(0,b,(nz))
    for i in n :
        for j in m:
            for k in q:
                Ux=(-P*j/(6*E*I))*(((6*l-3*i)*i)+(2+v)*(j**2+(h**2)/4))
                Uy=(P/(6*E*I))*(((3*v*(j**2))*(l-i))+((4+5*v)*((h**2)*i/4))+(3*l-i)*(i**2))
                Uz=0
                disp.append([i,j,k,Ux,Uy,Uz])
    return disp
