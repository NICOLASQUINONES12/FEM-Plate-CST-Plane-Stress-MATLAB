%We create the local stiffness matrix
%Counterclockwise system 

function local_stiffness_matrix = local_stiffness_matrix(x0,y0,x1,y1,x2,y2,E,v,t) 
    B = strain_displacement_matrix(x0,y0,x1,y1,x2,y2);
    D = constitutive_matrix(E,v);
    A = two_area_triangle(x0,y0,x1,y1,x2,y2)/2;
    k = t*A*(B'*D*B);
    local_stiffness_matrix = k;
end

function strain_displacement_matrix = strain_displacement_matrix(x0,y0,x1,y1,x2,y2)
    bheta_0 = y1-y2;
    bheta_1 = y2-y0;
    bheta_2 = y0-y1;
    theta_0 = x2-x1;
    theta_1 = x0-x2;
    theta_2 = x1-x0;
    two_Area = two_area_triangle(x0,y0,x1,y1,x2,y2);
    B = zeros(3,6);
    B(1,1) = bheta_0; B(1,3) = bheta_1; B(1,5) = bheta_2;
    B(2,2) = theta_0; B(2,4) = theta_1; B(2,6) = theta_2;
    B(3,1) = theta_0; B(3,2) = bheta_0; B(3,3) = theta_1; B(3,4) = bheta_1; B(3,5) = theta_2; B(3,6) = bheta_2;
    B = 1/(two_Area)*B;
    strain_displacement_matrix = B; 
end

%Constitutive_matrix for plane stresses
function constitutive_matrix = constitutive_matrix (E,v) 
    D = zeros(3,3);
    D(1,1) = 1; D(1,2) = v;
    D(2,1) = v; D(2,2) = 1;
    D(3,3) = (1-v)/2;
    D = E/(1-v^2)*D;
    constitutive_matrix = D;
end

%With a Counterclockwise system the area is positive
function two_area_triangle = two_area_triangle(x0,y0,x1,y1,x2,y2)
    two_area_triangle = x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1);
end