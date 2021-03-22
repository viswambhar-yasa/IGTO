#AUTHOR : YASA VISWAMBHAR REDDY
#MATRICULATION NUMBER : 65074
#Personal Programming Project
#--------------------------------------------------------------------------------------------#
#ANALYTIAL SOLUTION - Python file used to generate analytical solution for cantilever beam 2D
#--------------------------------------------------------------------------------------------#

from numpy import linspace,array

def exact_displacements(P,l,b,h,E,v,nx,ny,nz):
    '''
    function gives analytical solution for a 2D cantilever beam of length l and height h with point load acting at the top corner at the free end.

    Parameters
    ----------
    P : int 
        load or point force acting at the free end.

    l,b,h :int
        lenght, breath, width of the cantilever beam

    E :int
        Young's modulus

    v: float
        poission's ratio

    nx,ny,nz : int
            Number of element in x,y and z direction.
    '''
    #Moment of inertia is calculated
    I=b*(h**3)/12
    # to store the analytical displacement
    disp=[]
    #initialization of the mesh
    n=linspace(0,l,(nx))
    m=linspace(0,h,(ny))
    q=linspace(0,b,(nz))
    #looped over number mesh points
    for i in n :
        for j in m:
            for k in q:
                Ux=(-P*j/(6*E*I))*(((6*l-3*i)*i)+(2+v)*(j**2+(h**2)/4))                      #based on equ. 6.8
                Uy=(P/(6*E*I))*(((3*v*(j**2))*(l-i))+((4+5*v)*((h**2)*i/4))+(3*l-i)*(i**2))  #based on qua. 6.9
                Uz=0
                disp.append([i,j,k,Ux,Uy,Uz]) # displacement and their co-ordinates are stored for evaluation with numerical results.
    return disp
