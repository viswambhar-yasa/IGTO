from numpy import linspace,array


def exact_displacements(P,l,b,h,E,v,nx,ny):
    I=b*(h**3)/12
    disp=[]
    n=linspace(0,l,(nx+1))
    m=linspace(0,h,(ny+1))
    for i in n :
        for j in m:
            Ux=(-P*j/(6*E*I))*(((6*l-3*i)*i)+(2+v)*(j**2+(h**2)/4))
            Uy=(-P/(6*E*I))*(((3*v*(j**2))*(l-i))+((4+5*v)*((h**2)*i/4))+(3*l-i)*(i**2))
            disp.append([Ux,Uy])
    return disp
'''
length=48
height=12
width=1
option=0
nx=20
ny=6
nz=3
density=7850
volume_frac=0.5
pmax=3
rmin=1.5
load=-10

XI_DEGREE=2
ETA_DEGREE=2
NETA_DEGREE=2
    
Youngs_modulus=100000
poission_ratio=0.3

ex_dis=array(exact_displacements(P,l,b,h,E,v,nx,ny))
print(ex_dis[:,1:])

I=(w*(h**3))/12
disp=(P*(l**3))/(3*E*I)+((P*l)/((E/2*(1-v))*(5/6)*w*h))
print(disp)

i=1
j=1
I=w*(h**3)/12

Ux=(-P*j/(6*E*I))*(((6*l-3*i)*i)+(2+v)*(j**2+(h**2)/4))
Uy=(-P/(6*E*I))*(((3*v*(j**2))*(l-i))+((4+5*v)*((h**2)*i/4))+(3*l-i)*(i**2))

print(Uy,Ux)
'''
'''
I=h**3/12
x=l
y=h
Ux=-(P/(3*E*I))*((6*l-3*x)*x+(2+v)*((y**2)-((h**2)/4)))
Uy=(P/(6*E*I))*((3*v*(y**2)*(l-x))+((4+5*v)*((h**2)*(x/4)))+((3*l-x)*(x**2)))
print(Ux)
print(Uy)
'''