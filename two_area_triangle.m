%With a Counterclockwise system the area is positive
function two_area = two_area_triangle(x0,y0,x1,y1,x2,y2)
    two_area = x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1);
end