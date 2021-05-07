# Isogeometric analysis based topology optimization (IGTO)
<p align="center">
<img src="https://github.com/viswambhar-yasa/IGTO/blob/main/ezgifcom-gif-maker.gif"  width="500" height="350"/>
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

<p align="center">
<img src="https://github.com/viswambhar-yasa/Non-Linear-FEM/raw/master/program_structure.png"  />
</p>



#### Material Routine 
The material is consider to be visco-elastic and thus non-linearity is implemented by using Newton Raphson method.

#### Program Structure

<p align="center">
<img src="https://github.com/viswambhar-yasa/Non-Linear-FEM/raw/master/program_structure.png"  />
</p>


This is a repository of the files belonging to the assignment as a part of the module "Nonlinear Finite Element Methods" in the curriculum of Computational Materials Science.
