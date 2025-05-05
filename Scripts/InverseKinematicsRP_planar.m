function [q1,q2] = InverseKinematicsRP_planar(x,y,L)
%INVERSEKINEMATICSRP_PLANAR Summary of this function goes here
%   Detailed explanation goes here

    q2 = sqrt(x^2+y^2)-L;
    q1 = atan2(y,x);
fprintf("q1: %f, q2: %f\n", q1, q2);
end

