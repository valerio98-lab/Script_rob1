function [axis, theta] = RodriguezInverseSolution(R)

tol = 1e-6;

if (abs(det(R) - 1) > tol) || (norm(R' * R - eye(3)) > tol)
    error('Is not a valid rotation matrix. Watch Out!');
end

t = trace(R);
cosTheta = (t - 1) / 2;
cosTheta = max(min(cosTheta, 1), -1);

sinTheta = norm([R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)]) / 2;
fprintf('Il seno: %d \n', sinTheta)
theta = atan2(sinTheta, cosTheta);

if sinTheta > tol
    axis = 1/(2*sinTheta) * [R(3,2)-R(2,3);
                             R(1,3)-R(3,1);
                             R(2,1)-R(1,2)];
else
    if abs(theta) < tol
        disp('There is no axis-solution for theta=0. The axis returned is randomly chosen')
        axis = [1; 0; 0];
    else
        x = sqrt(max((R(1,1) + 1) / 2, 0));
        y = sqrt(max((R(2,2) + 1) / 2, 0));
        z = sqrt(max((R(3,3) + 1) / 2, 0));
        
        if abs(x) >= tol
            if R(1,2) < 0
                y = -y;
            end
            if R(1,3) < 0
                z = -z;
            end
        elseif abs(y) >= tol
            if R(1,2) < 0
                x = -x;
            end
            if R(2,3) < 0
                z = -z;
            end
        else
            if R(1,3) < 0
                x = -x;
            end
            if R(2,3) < 0
                y = -y;
            end
        end
        
        axis = [x; y; z];
        axis = axis / norm(axis);
    end
end

end

