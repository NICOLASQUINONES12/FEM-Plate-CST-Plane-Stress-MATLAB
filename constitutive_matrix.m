%Constitutive_matrix for plane stresses
function C = constitutive_matrix (E,v) 
    C = zeros(3,3);
    C(1,1) = 1; C(1,2) = v;
    C(2,1) = v; C(2,2) = 1;
    C(3,3) = (1-v)/2;
    C = E/(1-v^2)*C;
end
