%We create the local stiffness matrix
%Counterclockwise system 

function k = local_stiffness_matrix(x0,y0,x1,y1,x2,y2,E,v,t) 
    B = strain_displacement_matrix(x0,y0,x1,y1,x2,y2);
    C = constitutive_matrix(E,v);
    A = two_area_triangle(x0,y0,x1,y1,x2,y2)/2;
    k = t*A*(B'*C*B);
end
