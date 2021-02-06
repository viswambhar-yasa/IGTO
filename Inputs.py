TGREEN = '\033[32;1m'
TYELLOW =  '\033[33;1m' 
ENDC = '\033[m'
print(TGREEN+'\n Press Enter for default values or Enter values accordingly \n'+ENDC)
print('Ft_op=0 # OC, Ft_op=1 #MMA \n')
print('BC_op=0 #c-c-c-c, BC_op=1 #s-s-s-s, BC_op=2 #c-c-c-c load at end ,BC_op=3 #c-c-c-c analytical ,BC_p=4 #c-c-c-c two forces \n')
print(TYELLOW+ 'l h w nx ny nz load volume_fra penal rmin E v density BC_op Ft_op verbose' +ENDC)
Inputs = input()
Input_list=Inputs.split()
if len(Input_list)==0:
    length=12
    height=4
    width=3
    option=3
    nx=3
    ny=3
    nz=3
    density=7850
    volume_frac=0.5
    penal=3.5
    rmin=1.5
    load=-200
    mesh_disp=False
    iterative_display=False
    optimizer='MMA'
    Youngs_modulus=100000
    poission_ratio=0.35

elif len(Input_list)==16:    
    #print(inp)
    length=float(Input_list[0])
    height=float(Input_list[1])
    width=float(Input_list[2]) 
    option=int(Input_list[13]) 
    nx=int(Input_list[3]) 
    ny=int(Input_list[4]) 
    nz=int(Input_list[5]) 
    volume_frac=float(Input_list[7]) 
    penal=float(Input_list[8]) 
    rmin=float(Input_list[9]) 
    load=float(Input_list[6]) 
    density=float(Input_list[12]) 
    verbose=int(Input_list[15])  
    optimizer=int(Input_list[14]) 
    Youngs_modulus=float(Input_list[10])
    poission_ratio=float(Input_list[11])
    if optimizer==0:
        optimizer='OC'
    else:
        optimizer='MMA'

    if verbose==0:
        mesh_disp=False
        iterative_display=False
    else:
        mesh_disp=False
        iterative_display=True
else: 
    raise ValueError('No enough aruguments')
#print(l h w nx ny nz load volume_fra penal rmin E v density BC_op Ft_op verbose)
#48 12 1 13 11 3 -200 0.5 3.5 1.25 100000 0.3 7850 3 1 0 
#12 3 4 13 3 9 -200 0.15 10 5 150000 0.35 7850 3 1 1 
