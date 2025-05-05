function R_finale = RodriguezMatrix(r, theta)
    r = r';
    r = r/norm(r);
    
    skew_sym = @(r) [0, -r(3), r(2);
                     r(3), 0, -r(1); 
                    -r(2), r(1), 0];
    
    Rodriguez = @(r, theta) r*r' + (eye(3) - r*r')*cos(theta) + skew_sym(r)*sin(theta);
    R_finale = Rodriguez(r, theta);

    
    
end

