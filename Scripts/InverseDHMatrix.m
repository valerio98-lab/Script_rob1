function [alpha, theta, a ,d] = InverseDHMatrix(R, verbose);
%INVERSEDHMATRIX Return the values of DH-convention w.r.t the matrix R
%  

alpha = atan2(R(3,2), R(3,3));
theta = atan2(R(1,2), R(1,1));
a = R(1,4) * cos(theta) + R(2,4)*sin(theta);
d = R(3,4);

if verbose
    fprintf('alpha %f \n', alpha)
    fprintf('theta %f \n', theta)
    fprintf('a %f \n', a)
    fprintf('d %f \n', d)

end

