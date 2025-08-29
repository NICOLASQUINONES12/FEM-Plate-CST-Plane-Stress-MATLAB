% ---
% jupyter:
%   jupytext:
%     formats: ipynb,m:percent
%     text_representation:
%       extension: .m
%       format_name: percent
%       format_version: '1.3'
%       jupytext_version: 1.17.2
%   kernelspec:
%     display_name: Matlab
%     language: matlab
%     name: matlab
% ---

% %% [markdown]
% # A Finite Element Solution for a plate with MATLAB code asumming plane stress. 

% %% [markdown]
% The problem is show in the figure 6-16. 

% %% [markdown]
% ![Figure 1](figure1.JPG)

% %% [markdown]
% This problem comes from the chapter 6.5 of the book: *A firts Course in the Finite Element Method* written by Daryl L. Logan. In this 
% chapter a solution is make by subdivided the plate in two CST. The purpose of this project is to make a code that can performed more divisions. 

% %% [markdown]
% ## Input Variables

% %% [markdown]
% Here you can change the variables. 

% %%
% Geometry definition 

%The mesh is divided in n_x*n_y squares, n_x division in the x direction and n_y divisions in the y direction. 
%At each rectangle the CST is generated joining the lower-left node with the upper-right node of the rectangle. 

t = 1; % in
E = 30*10^6; %psi
v = 0.3; 
n_x = 10; 
n_y = 10;
%For each element there are additional divisions to calculate the internal stresses and strains
n_element_x = 5; 
n_element_y = 5;
x0 = 0;
y0 = 0;
x1 = 20;
y1 = 10; 
%Surface Forces
p = 1000; %Asumming that positive values are for traction and negative values for comprension 
%At the end of the code is observed a graphic of deformation of the nodes, with the following variable you can change the displacement scale in order 
% to see the difference between the original nodes and the displaced nodes more clearly. 
displacement_scale = 1000;


% %% [markdown]
% ## Important notes of the FEM with CST elements

% %% [markdown]
% - When a large number of CST elements is considered, the model becomes stiffer than it should be. At the following section, the displacement in $x$ and the strain in $x$ as a single element model are shown. It can be observed that for large $n$ values such as 20 or 25 CST elements, the displacement in $x$ and the strain in $x$ are lower than the solution with the single-element model.
% - It is recommended to use 10 divisions in x and y to have an accurate solution.

% %% [markdown]
% ## Solution as approximation of one element 

% %% [markdown]
% If one element is used, the displacement is going to be:
% $$
% \delta = \frac{P L}{E A}
% $$
% This result can be used to compare the values of the CST method with one element method. This value is situable for displacements in $x$ ($U_{xy}$) in the graphics. 

% %%
A = (y1-y0)*t;
L = (x1-x0);
displacement_x_one_element = p*A*L/(E*A);
displacement_x_one_element

% %% [markdown]
% The strain is equal to: 
% $$
% \epsilon = \frac{\delta}{L}
% $$
% This value is situable for strain in x ($\epsilon_{x}$) in the graphics.  

% %%
strain_x_one_element = displacement_x_one_element/L; 
strain_x_one_element

% %% [markdown]
% ## Solution as approximation of FEM with CST elements 

% %% [markdown]
% The following approach to solve the problem is with Constant Strain Triangles (CST). Plane stress were assumed so the stresses in the z direction are cero.
%
% The following image shows the CST mesh for this assembly.

% %%
%This code is to divided the plate in x and y. It must be performed once and it is valid for all the code. 
x = linspace(x0,x1,n_x+1); 
y = linspace(y0,y1,n_y+1); 
%-------------------------------------------------%
if (n_x<=25) || (n_y<=25)
    [X_gr, Y_gr] = meshgrid(x,y);
    % Horizontal lines 
    for i = 1:size(X_gr,1)
        plot(X_gr(i,:),Y_gr(i,:),'k');
        hold on
    end
    % Vertical lines 
    for j = 1:size(X_gr,2)
        plot(X_gr(:,j),Y_gr(:,j),'k');
        hold on
    end
    % Diagonal lines 
    for i = 1:(size(X_gr,1)-1)
        for j = 1:(size(X_gr,2)-1)
            x0 = X_gr(i,j);
            x1 = X_gr(i+1,j+1);
            y0 = Y_gr(i,j);
            y1 = Y_gr(i+1,j+1);
            x_point = [x0,x1];
            y_point = [y0,y1];
            plot(x_point, y_point, 'k');
            hold on
        end
    end
    title('CST elements');
    xlabel('x coordinates');
    ylabel('y coordinates');
    ax = gca;
    ax.TickLength = [0 0];
else
    disp('The image of CST elements can be shown only if n_x <= 25 and n_y <= 25.');
    disp('n_x are number of divisions in x and n_y are number or divisions in y.');
end

% %% [markdown]
% ## Formulation, Assembly and Solution

% %%
%Numeration of nodes
%This function enumerate the nodes from left to right starting in the coordinate (1,1)
node_numeration = @(i,j) (j-1)*(n_x+1)+i; 

%Global matrix
n_nodes = (n_x+1)*(n_y+1); 
ndfo = 2*n_nodes; 
K = sparse(ndfo,ndfo);

%Vector of global displacements 
D = NaN(ndfo,1);
%Vector of global Forces; 
Q = NaN(ndfo,1);

%Surface force 
%The surface force is divided in each element so the force in each element is (p*L*t), then it is transmitted 
%to the nodes with the value (pLt/2) to each node. 

%Value of L
%Because of the problem L es divided in the y direction 
L = (y1-y0)/n_y;
%Node force 
node_force = p*L*t/2;

%Application of the boundary conditions and assembly of the global matrices Q, D and K. 

for i = 1:n_x
    for j = 1:n_y
    
        %Application of the boundary conditions 
        %--------------------------------------------------------------------------------%
        %Displacements are cero in the nodes (1,1) to (1,j+1) - left part of the body. 
        if i == 1
            n0 = node_numeration(i,j);
            dfo0 = [2*n0-1,2*n0];
            d = [0;0];
            D(dfo0) = d;
            %Additionally:
            if j == n_y
               %node (1,j+1) 
               n1 = node_numeration(i,j+1);
               dfo1 = [2*n1-1,2*n1]; 
               D(dfo1) = d;
            end
        end 
        
        %Forces in the right side of the plate 
        %Calculation of surface forces (from nodes (i+1,1) to (i+1,j+1)) and assembly to the vector of global Forces 
        if i == n_x
            n0 = node_numeration(i+1,j);
            dfo0 = [2*n0-1,2*n0];
            q = [node_force;0];
            
            if j == 1
                Q(dfo0) = q;  
            else 
                Q(dfo0) = 2*q;
            end
            
            %Additionally:  
            if j == n_y
                %node (i+1,j+1)
                n1 = node_numeration(i+1,j+1);
                dfo1 = [2*n1-1,2*n1]; 
                Q(dfo1) = q;
            end    
        end
        
        %Forces are cero from the nodes (2,1) to (i,j). 
        %Cero forces
        if (i>1)&&(i<=n_x)
            n0 = node_numeration(i,j);
            df0 = [2*n0-1,2*n0];
            q = [0;0];
            Q(df0) = q;
        end
        
        %Additionally
        if (i>1)&&(j == n_y)
            n0 = node_numeration(i,j+1);
            df0 = [2*n0-1,2*n0];
            q = [0;0];
            Q(df0) = q;
        end
        
        %The values that appear with NaN in the matrices Q and D are incognites, the values that are already known were replaced in the code above. 
        
        %Calculation of local stiffness matrix and assembly to global stiffness matrix 
        %--------------------------------------------------------------------------------%
        
        %The local stiffness matrix was coded in the file: "local_stiffness_matrix.m"
        
        %Local nodes transformed to global nodes numeration 
        n1 = node_numeration(i,j);
        n2 = node_numeration(i+1,j);
        n3 = node_numeration(i+1,j+1);
        n4 = node_numeration(i,j+1); 
        
        %Nodes coordinates
        x1 = x(i); y1 = y(j);
        x2 = x(i+1); y2 = y(j);
        x3 = x(i+1); y3 = y(j+1);
        x4 = x(i); y4 = y(j+1);
        
        %Triangle1: n1(x1,y1), n3(x3,y3), n4(x4,y4) 
        k1 = local_stiffness_matrix(x1,y1,x3,y3,x4,y4,E,v,t);
        
        %Degrees of freedom numeration 
        dfo1 = [2*n1-1,2*n1,2*n3-1,2*n3,2*n4-1,2*n4];
        
        %Assembly
        K(dfo1,dfo1) = K(dfo1,dfo1) + k1;

        %Triangle2: n1(x1,y1), n2(x2,y2), n3(x3,y3) 
        k2 = local_stiffness_matrix(x1,y1,x2,y2,x3,y3,E,v,t);
        
        %Degrees of freedom numeration 
        dfo2 = [2*n1-1,2*n1,2*n2-1,2*n2,2*n3-1,2*n3];
        
        %Assembly
        K(dfo2,dfo2) = K(dfo2,dfo2) + k2;
    end
end

%Solution of the equation Q = K*D
%--------------------------------------------------------------------------------%

%The matrices are subdivided
%u is used for unknown and k for known 
D_k_positions = find(~isnan(D)); %The displacements that are imposed
D_u_positions = find(isnan(D)); %The displacements that are not known

K_uu = K(D_u_positions,D_u_positions);
K_uk = K(D_u_positions,D_k_positions);
K_kk = K(D_k_positions,D_k_positions);
K_ku = K(D_k_positions,D_u_positions);

D_u = D(D_u_positions);
D_k = D(D_k_positions);
Q_k = Q(D_u_positions);%Values of Q in the dfo where the displacements are not known 

D_u = K_uu\(Q_k-K_uk*D_k);
Q_u = K_kk*D_k+K_ku*D_u;

D(D_u_positions) = D_u;
Q(D_k_positions) = Q_u;

% %% [markdown]
% ## Graphics 

% %%
%In this code is calculate the displacements, strains and stresses inside each element
%--------------------------------------------------------------------------------%

%Global matrices 
U_x_y = zeros(n_y*n_element_y + 1, n_x*n_element_x + 1);
V_x_y = zeros(n_y*n_element_y + 1, n_x*n_element_x + 1);
strain_x = zeros(n_y*n_element_y + 1,n_x*n_element_x + 1);
strain_y = zeros(n_y*n_element_y + 1,n_x*n_element_x + 1);
strain_shear = zeros(n_y*n_element_y + 1,n_x*n_element_x + 1);
stress_x = zeros(n_y*n_element_y + 1,n_x*n_element_x + 1);
stress_y = zeros(n_y*n_element_y + 1,n_x*n_element_x + 1);
tau_x_y = zeros(n_y*n_element_y + 1,n_x*n_element_x + 1); 
%It is calculated only once: 
C = constitutive_matrix(E,v);

for i = 1:(n_x)
    for j = 1:(n_y)
        %The rectangle (i,j,i+1,j+1) is divided into two triangles
        %counterclockwise the area is positive
        %Node numeration in the small rectangle 
        node_0 = node_numeration(i,j);
        node_1 = node_numeration(i+1,j);
        node_2 = node_numeration(i,j+1);
        node_3 = node_numeration(i+1,j+1);
        
        %First triangle
        x0 = [x(i),x(i+1),x(i+1)];
        y0 = [y(j),y(j),y(j+1)];
        u0 = [D(2*node_0-1),D(2*node_1-1),D(2*node_3-1)];
        v0 = [D(2*node_0),D(2*node_1),D(2*node_3)];
        d0 = [u0(1);v0(1);u0(2);v0(2);u0(3);v0(3)];
        
        %Second triangle
        x1 = [x(i),x(i+1),x(i)];
        y1 = [y(j),y(j+1),y(j+1)];
        u1 = [D(2*node_0-1),D(2*node_3-1),D(2*node_2-1)];
        v1 = [D(2*node_0),D(2*node_3),D(2*node_2)];
        d1= [u1(1);v1(1);u1(2);v1(2);u1(3);v1(3)];
        
        %The limits for the rectangle
        xmin = min([x0,x1]);
        ymin = min([y0,y1]);
        xmax = max([x0,x1]);
        ymax = max([y0,y1]);
        local_x = linspace(xmin,xmax,n_element_x+1);
        local_y = linspace(ymin,ymax,n_element_y+1);
        
        %Local matrices
        %local displacements
        u_local = zeros(length(local_y),length(local_x));
        v_local = zeros(length(local_y),length(local_x));
        %local strains
        strain_x_local = zeros(length(local_y),length(local_x));
        strain_y_local = zeros(length(local_y),length(local_x));
        strain_shear_local = zeros(length(local_y),length(local_x));
        %local stresses
        stress_x_local = zeros(length(local_y),length(local_x));
        stress_y_local = zeros(length(local_y),length(local_x));
        tau_x_y_local = zeros(length(local_y),length(local_x));
        
        %First triangle
        B0 = strain_displacement_matrix(x0(1),y0(1),x0(2),y0(2),x0(3),y0(3));
        eps0 = B0*d0;
        sig0 = C*eps0;
        
        for k = 1:(n_element_y+1)   
            for m = 1:(n_element_x+1)
                %Verification that the point is inside the triangle
                lambda1 = barycentric_lambda(x0(1), y0(1), x0(2), y0(2), x0(3), y0(3), local_x(m), local_y(k));
                lambda2 = barycentric_lambda(x0(2), y0(2), x0(3), y0(3), x0(1), y0(1), local_x(m), local_y(k));
                lambda3 = barycentric_lambda(x0(3), y0(3), x0(1), y0(1), x0(2), y0(2), local_x(m), local_y(k));
                if (lambda1<0) || (lambda2<0) || (lambda3<0)
                    continue; 
                end
                
                %to calculate the local displacements inside the triangle. The shape functions had been used. 
                u_local(k,m) = u_x_y(u0(1),u0(2),u0(3),x0(1),y0(1),x0(2),y0(2),x0(3),y0(3),local_x(m), local_y(k));
                v_local(k,m) = u_x_y(v0(1),v0(2),v0(3),x0(1),y0(1),x0(2),y0(2),x0(3),y0(3),local_x(m), local_y(k));
                
                % to calculate the strains 
                strain_x_local(k,m) = eps0(1);
                strain_y_local(k,m) = eps0(2);
                strain_shear_local(k,m) = eps0(3);
                
                %to calculate the stresses
                stress_x_local(k,m) = sig0(1);
                stress_y_local(k,m) = sig0(2);
                tau_x_y_local(k,m) = sig0(3);
            end
        end
        
        %Second triangle
        B1 = strain_displacement_matrix(x1(1),y1(1),x1(2),y1(2),x1(3),y1(3));
        eps1 = B1*d1;
        sig1 = C*eps1;
        
        for k = 1:(n_element_y+1)     
            for m = 1:(n_element_x+1)
                %Verification that the point is inside the triangle
                lambda1 = barycentric_lambda(x1(1), y1(1), x1(2), y1(2), x1(3), y1(3), local_x(m), local_y(k));
                lambda2 = barycentric_lambda(x1(2), y1(2), x1(3), y1(3), x1(1), y1(1), local_x(m), local_y(k));
                lambda3 = barycentric_lambda(x1(3), y1(3), x1(1), y1(1), x1(2), y1(2), local_x(m), local_y(k));
                if (lambda1<0) || (lambda2<0) || (lambda3<0)
                    continue; 
                end
                
                %to calculate the local displacements inside the triangle. The shape functions had been used. 
                u_local(k,m) = u_x_y(u1(1),u1(2),u1(3),x1(1),y1(1),x1(2),y1(2),x1(3),y1(3),local_x(m), local_y(k));
                v_local(k,m) = u_x_y(v1(1),v1(2),v1(3),x1(1),y1(1),x1(2),y1(2),x1(3),y1(3),local_x(m), local_y(k));

                % to calculate the strains 
                strain_x_local(k,m) = eps1(1);
                strain_y_local(k,m) = eps1(2);
                strain_shear_local(k,m) = eps1(3);
                
                %to calculate the stresses
                stress_x_local(k,m) = sig1(1);
                stress_y_local(k,m) = sig1(2);
                tau_x_y_local(k,m) = sig1(3);
            end
        end
        
        %Assembly
        row_start = (j-1)*n_element_y + 1;
        row_end = j*n_element_y + 1;
        col_start = (i-1)*n_element_x + 1;
        col_end = i*n_element_x + 1;
        
        U_x_y(row_start:row_end, col_start:col_end) = u_local;
        V_x_y(row_start:row_end, col_start:col_end) = v_local;
        strain_x(row_start:row_end, col_start:col_end) = strain_x_local;
        strain_y(row_start:row_end, col_start:col_end) = strain_y_local;
        strain_shear(row_start:row_end, col_start:col_end) = strain_shear_local;
        stress_x(row_start:row_end, col_start:col_end) = stress_x_local;
        stress_y(row_start:row_end, col_start:col_end) = stress_y_local;
        tau_x_y(row_start:row_end, col_start:col_end) = tau_x_y_local;
    end
end

%Figures plot

x_plot = linspace(min(x),max(x),n_x*n_element_x+1);
y_plot = linspace(min(y),max(y),n_y*n_element_y+1);
[X,Y] = meshgrid(x_plot, y_plot);



% %%
%Importan annotation: When it is used pcolor to graph the values of the matrices U_x_y, V_x_y, etc., the values at the final row
%and final column are lost because the color is taken from the left lower side of the grid. 
%--------------------------------------------------------------------------------%
figure;
pcolor(X, Y, U_x_y);
axis equal tight;
shading flat; 
colorbar;
title('U_x_y Distribution');
xlabel('X coordinate');
ylabel('Y coordinate');

% %%
%Graphic V_x_y
figure;
pcolor(X, Y, V_x_y);
axis equal tight;
shading flat;
colorbar;
title('V_x_y Distribution');
xlabel('X coordinate');
ylabel('Y coordinate');

% %%
%Graphic strain_x
figure;
pcolor(X, Y, strain_x);
axis equal tight;
shading flat;
colorbar;
title('strain_x');
xlabel('X coordinate');
ylabel('Y coordinate');

% %%
%Graphic strain_y
figure;
pcolor(X, Y, strain_y);
axis equal tight;
shading flat;
colorbar;
title('strain_y');
xlabel('X coordinate');
ylabel('Y coordinate');

% %%
%Graphic shear_strain
figure;
pcolor(X, Y, strain_shear);
axis equal tight;
shading flat;
colorbar;
title('shear strain');
xlabel('X coordinate');
ylabel('Y coordinate');

% %%
%Graphic stress_x
figure;
pcolor(X, Y, stress_x);
axis equal tight;
shading flat;
colorbar;
title('stress_x');
xlabel('X coordinate');
ylabel('Y coordinate');

% %%
%Graphic stress_y
figure;
pcolor(X, Y, stress_y);
axis equal tight;
shading flat;
colorbar;
title('stress_y');
xlabel('X coordinate');
ylabel('Y coordinate');


% %%
%Graphic shear stress
figure;
pcolor(X, Y, tau_x_y);
axis equal tight;
shading flat;
colorbar;
title('shear stress');
xlabel('X coordinate');
ylabel('Y coordinate');

% %% [markdown]
% # Deformation of the plate 

% %%
%Here the deformation of the plate is observed through the displacement of the nodes. 
%Original Coordinates
[X_grid,Y_grid] = meshgrid(x,y);
X_grid_tra = X_grid';
Y_grid_tra = Y_grid';
x_points = X_grid_tra(:);
y_points = Y_grid_tra(:);
%Displacements of u for x with corresponds to the odd dfos and v for y with even dfos. 
u = D(1:2:end);
v = D(2:2:end);

%Displacements in x and y 
x_displaced = x_points +u*displacement_scale;
y_displaced = y_points +v*displacement_scale;

%Figures 
figure;
hold on; 
% original nodal points and displaced nodal points 
plot(X_grid(:), Y_grid(:), 'bo', 'MarkerFaceColor', 'b');
plot(x_displaced, y_displaced, 'ro', 'MarkerFaceColor', 'r');
hold off; 

%Titles, axes titles, legends and texts

h = title('Original Nodes vs Deformed Nodes', 'Units', 'Normalized');
set(h, 'Position', [0.5, 1.05, 0]);

xlabel('X coordinates','FontSize', 12);
ylabel('Y coordinates','FontSize', 12);
legend('Original nodes','Deformed nodes');
y_text_position = min(y_points) - (max(y_points) - min(y_points)) * 0.1; 
text(min(x_points+1), y_text_position, sprintf('Displacement Scale: %d', displacement_scale), 'FontSize', 10, 'Interpreter', 'none');

% %% [markdown]
% ## Matlab version 

% %%
version
