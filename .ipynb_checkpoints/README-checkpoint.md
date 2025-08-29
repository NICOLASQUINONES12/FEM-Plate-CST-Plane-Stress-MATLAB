# Implementation of FEM for a 2D Plate using CST Elements under the Plane Stress Assumption (MATLAB)

This project implements the Finite Element Method (FEM) for 2D plates using Constant Strain Triangle (CST) elements under the plane stress assumption.  
The solution corresponds to the problem illustrated in figure 6.16.
![Figure 1](figure1.JPG)

This problem comes from the chapter 6.5 of the book:  
*A First Course in the Finite Element Method* by **Daryl L. Logan**.

## Features

- FEM implementation with CST elements.  
- Visualization of displacements, stresses, and strains of the plate.  
- Scalable visualization of node deformations.  
- Comparison of the solution using a single-element model.

## Requirements

- MATLAB R2025a Update 1 or later

## Usage

Option 1: MATLAB
Open the `fem_cst.m` script directly in MATLAB and run it. This will execute the code within the MATLAB environment.

Option 2: Jupyter Notebook 
Open the `fem_cst.ipynb` file in a Jupyter Notebook. This requires that you have a MATLAB kernel configured in your Jupyter environment to execute the code.

## Example

For a 10×10 mesh with the following values:
- t = 1;            % thickness [in]
- E = 30*10^6;      % Young’s modulus [psi]
- v = 0.3;          % Poisson’s ratio
- n_x = 10;
- n_y = 10;
- n_element_x = 5;
- n_element_y = 5;
- x0 = 0;
- y0 = 0;
- x1 = 20;
- y1 = 10;
- p = -1000;        % Load [psi]
- displacement_scale = 1000;

Some Expected outputs include:
- Graphic of CST Elements in the plate
- Displacement in x (U_xy)
- Displacement in y (V_xy)
- Strain in x
- Strain in y
- Shear strain
- Stress in x
- Stress in y
- Shear stress 
- Deformed mesh grid 

Example results:

Displacements in x
![Displacement in x](figure1_u_x_y_example.JPG)
Strain in x
![Strain in x.](figure2_strain_x_example.JPG)
Deformed mesh grid 
![Deformed mesh grid.](figure_3_deformed_nodes_examples.JPG)

## Repository Structure

- `fem_cst.ipynb`: The main Jupyter Notebook for the project, demonstrating the FEM analysis.
- `fem_cst.m`: An equivalent MATLAB script of the main FEM program.
- `barycentric_lambda.m`: A MATLAB function to compute barycentric coordinates for a triangle.
- `constitutive_matrix.m`: A MATLAB function for generating the material's constitutive matrix.
- `local_stiffness_matrix.m`: A MATLAB function that calculates the local stiffness matrix.
- `strain_displacement_matrix.m`: A MATLAB function that defines the strain-displacement matrix B.
- `two_area_triangle.m`: A MATLAB utility function to compute the area of a CST element. 
- `u_x_y.m`: A MATLAB function that determines displacements inside a given element.
- `figure1.JPG`: An image illustrating the initial problem setup.
- `figure1_u_x_y_example.JPG`: An example image showing the displacement results.
- `figure2_strain_x_example.JPG`: An example image showing the strain results.
- `figure3_deformed_nodes_example.JPG`: An example image of the deformed mesh nodes.

## Reproducibility
To reproduce the results of this project, follow these steps:

1. Clone this repository: Use Git to get a local copy of all the project files.
2. Open the MATLAB script or notebook: Open fem_cst.m directly in MATLAB or fem_cst.ipynb in a Jupyter environment with a MATLAB kernel.
3. Run all cells: Execute the entire script or notebook sequentially.
4. Compare results: The generated output (figures, values) should match the reference images included in the repository.
  
## Author
Developed by Nicolas Quiñones 