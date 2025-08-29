%If barycentric_lambda>=0 the point (x,y) is inside the triangle
function lambda = barycentric_lambda(x0,y0,x1,y1,x2,y2,x,y)
    % Calculates the lambda associate with the vertice P(x0,y0)
    two_area = two_area_triangle(x0,y0,x1,y1,x2,y2);
    lambda = ((y1 - y2)*(x - x2) + (x2 - x1)*(y - y2)) / two_area;
end
