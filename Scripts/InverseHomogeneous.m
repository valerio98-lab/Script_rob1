function R = InverseHomogeneous(T)
    R_rot = T(1:3,1:3);
    p = T(1:3, 4);

    R = [R_rot', -R_rot'*p; 0, 0, 0, 1];
end

