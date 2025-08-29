# Implementation of FEM for a 2D Plate using CST Elements under the Plane Stress Assumption (MATLAB)

This project implements the Finite Element Method (FEM) for 2D plates using Constant Strain Triangle (CST) elements under the plane stress assumption.  
The solution corresponds to the problem illustrated in figure 6.16.
![Figure 1](figure1.JPG)

This problem is based on chapter 6.5 of the book:  
*A First Course in the Finite Element Method* by **Daryl L. Logan**.

## Features

- FEM implementation with CST elements.  
- Visualization of displacements, stresses, and strains of the plate.  
- Scalable visualization of node deformations.  
- Comparison of the solution using a single-element model.

## Requirements

- MATLAB R2025a Update 1 or later

## Usage

1. Open MATLAB.  
2. Open the script `fem_cst.m` or the Jupyter Notebook `fem_cst.ipynb` using MATLAB.  
3. Run all cells to reproduce the solution.

## Example

For a 10×10 mesh with the following values:
t = 1;            % thickness [in]
E = 30*10^6;      % Young’s modulus [psi]
v = 0.3;          % Poisson’s ratio
n_x = 10; 
n_y = 10;
n_element_x = 5; 
n_element_y = 5;
x0 = 0; 
y0 = 0;
x1 = 20; 
y1 = 10;
p = -1000;        % Load [psi]
displacement_scale = 1000;

Expected outputs include:

- Displacement in x.
- Strain in x.
- Deformed mesh grid. 

Example results:
![Displacement in x](figure1_u_x_y_example.JPG)
![Strain in x.](figure2_strain_x_example.JPG)
![Deformed mesh grid.](figure_3_deformed_nodes_examples.JPG)

Repository Structure

├── fem_cst.ipynb                # Main Jupyter Notebook
├── fem_cst.m                    # Equivalent MATLAB script
├── barycentric_lambda.m         # Function for barycentric lambda for a triangle (MATLAB)
├── constitutive_matrix.m        # Function for constitutive matrix (MATLAB)
├── local_stiffness_matrix.m     # Function for local stiffness matrix (MATLAB)
├── strain_displacement_matrix.m # Function for strain-displacement matrix B (MATLAB)
├── two_area_triangle.m          # Function to compute area (MATLAB)
├── u_x_y.m                      # Displacements inside a element, returns a value (MATLAB)
├── figure1.JPG                  # Example problem image
├── figure1_u_x_y_example.JPG    # Example result (displacement)
├── figure2_strain_x_example.JPG # Example result (strain)
├── figure3_deformed_nodes_example.JPG # Example result (deformed nodes)
└── README.md                     # Project documentation

