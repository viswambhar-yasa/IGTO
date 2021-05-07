<p align="center">

# IGTO
Isogeometric analysis based topology optimization
<\p>
The main idea behind structural optimization is to get structures with minimum
weight for different stresses and design constraints. Topology optimization is a
mathematical method of optimizing the distribution of material within the domain
so that the structure would satisfy the design parameters. The design parameters
include geometry of the domain, applied load and position of the load, the
amount of mass which has remain after optimization. We implemented SIMP(Solid
Isotropic Material with Penalization) method in order to remove porosity as this
method is quite effective for minimum weight topology

### Problem Statement 
The problem of creep of a thick-walled pipe under internal pressure p is considered
as sketched . The pressure rises linearly up to its final value pmax and
is then hold until final time. Plain strain along z direction is 0 conditions are
assumed.

<p align="center">
<img src="https://github.com/viswambhar-yasa/Non-Linear-FEM/raw/master/Annotation 2020-09-06 174514.png"  />
</p>

#### Material Routine 
The material is consider to be visco-elastic and thus non-linearity is implemented by using Newton Raphson method.

#### Program Structure

<p align="center">
<img src="https://github.com/viswambhar-yasa/Non-Linear-FEM/raw/master/program_structure.png"  />
</p>


This is a repository of the files belonging to the assignment as a part of the module "Nonlinear Finite Element Methods" in the curriculum of Computational Materials Science.
