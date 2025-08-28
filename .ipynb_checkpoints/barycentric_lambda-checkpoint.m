%lambda baricentrico
function lambda = barycentric_lambda(x0,y0,x1,y1,x2,y2,x,y)
    % Calculates the lambda associate with the vertice P(x0,y0)
    two_area = two_area_triangle(x0,y0,x1,y1,x2,y2);
    lambda = ((y1 - y2)*(x - x2) + (x2 - x1)*(y - y2)) / two_area;
end

%With a Counterclockwise system the area is positive
function two_area_triangle = two_area_triangle(x0,y0,x1,y1,x2,y2)
    two_area_triangle = x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1);
end