% u_x_y is how the displacements varies inside the element, u could be v, refering to x displacements and y displacements, respectively
%For this porpuse it is used the shape functions for a triangle shape. 

function displacement_x_y = u_x_y(u0,u1,u2,x0,y0,x1,y1,x2,y2,x,y)
    alpha_0 = x1 * y2 - y1 * x2;
    alpha_1 = y0 * x2 - x0 * y2;
    alpha_2 = x0 * y1 - y0 * x1;
    betha_0 = y1 - y2;
    betha_1 = y2 - y0;
    betha_2 = y0 - y1;
    gamma_0 = x2 - x1;
    gamma_1 = x0 - x2;
    gamma_2 = x1 - x0;    
    two_area = two_area_triangle(x0,y0,x1,y1,x2,y2);
    N_0 = 1/two_area*(alpha_0+betha_0*x+gamma_0*y);
    N_1 = 1/two_area*(alpha_1+betha_1*x+gamma_1*y);
    N_2 = 1/two_area*(alpha_2+betha_2*x+gamma_2*y);
    displacement_x_y = N_0*u0+N_1*u1+N_2*u2;
end
