# Isogeometric analysis based topology optimization (IGTO)
<p align="center">
<img src="https://github.com/viswambhar-yasa/IGTO/blob/main/ezgifcom-gif-maker.gif"/>
</p>

## INTRODUCTION

Topology optimization can be described as binary compliance problem which is
use to find a "black and white" layout that minimizes the compliance subjected to
volume constrain. There are many frameworks like homogenization method and
density-based approach.In our case, the material properties are assumed constant
within each element. We implement density-based approach as it restricts the
formation of pores and solved the minimum compliance effectively.
In density-based approach, the problem is parametrized by the material density
distribution. The stiffness tensor are determined using a power-law interpolation
function. The power-law implicitly penalizes the intermediate density values to
remove pores in structure layout. This penalization method is referred as Solid
Isotropic Material with Penalization (SIMP).

## USER MANUAL
The following steps should be performed to run the program and test cases. All
the files are written in python.
<p align="center">
<img src="https://github.com/viswambhar-yasa/IGTO/blob/main/Document/list_of_files.png" />
</p>

### Procedure to run the program
All python files have to be placed in same folder and working directory has to be
same to run the python program. Python version 3.7.4 can be used to run this
python code.
1. A file name main program.py is the starting point of the IGTO program.
command: `python main program.py`
(a) Required python external libraries like 'numpy', 'matplotlib', 'pyvista',
'pyEVTK','pytest' are checked and installed.
<p align="center">
<img src="https://github.com/viswambhar-yasa/IGTO/blob/main/Document/installed_lib.png" />
</p>
2. A prompt appears will asks the user to choose if time analysis has to be
performed or not.

0- Time analysis (log files are generated)

Enter - For normal execution without any log files.

3.Then the inputs have to be given or press enter to run default values.

l h w nx ny nz load volume_fra penal rmin E v density BC_op Ft_op
verbose

8 5 1 35 25 3 -100 0.4 3 1.5 150000 0.3 7850 0 1 1

### INPUTS WITH DEFAULT VALUES:

Length(l) : 8

Height(h) : 5

Width(w) : 1

nx(no of elements along x) : 35

ny(no of elements along y) : 25

nz(no of elements along z) : 3

load : -100

volume_fraction : 0.4

penal : 3

rmin : 1.5

E(Youngs modulus) : 150000

v(poisson ratio) : 0.3

density : 7850

BC_op : 0

Ft_op : 1

verbose : 1

#### Boundary condition option (BC op):
Options :

0- Cantilever beam with load along the bottom edge of the free end.

1- Simple supported with point load at bottom center.

2- Cantilever beam with point load at bottom of the free end.

3- Cantilever beam with point load at the free end (2d case
loading at y=height and x=length).

4- Cantilever beam with two forces at top and bottom of the free
end .

5- Cantilever beam with point load at the center of the free end.

6- Simple supported with load at top center.

#### Optimizer option (Ft op):
Options :

0-OC (optimality criterion)

1-MMA (method of moving asymptotes)

#### Verbose:
Options :

0-Will not print plots using pyvista only VTK file is generated.

1- Plots are generated and stored in results folder.

## Procedure to run the test case
We used pytest to perform unit, functional and patch testing. All test files are
the python files which contain respective test cases.
All files should be placed in the same folder and the working directory has to be
same, so that test cases can run.
<p align="center">
<img src="https://github.com/viswambhar-yasa/IGTO/blob/main/Document/list_of_tests.png"  />
</p>
1. test Inputs.py : Test cases for input parameters

2. test geometry.py : Test cases on shape function and assembly

3. test element routine.py : Test cases on integration scheme and sanity checks
on other element.

4. test optimization.py : Test cases on optimizer function OC(optimality cri-
terion) and MMA(method of moving asymptotes)

5. test rigid body motion.py : Test cases to validate global stiffness matrix
and boundary conditions by performing rigid body translation and rotation.

6. test patch.py : Test cases on constant stress patch test and comparing an-
alytical solution with numerical for lower and higher order shape functions.

To run all the test cases, copy all test files in same folder and enter `PYTEST`
command on the terminal.
