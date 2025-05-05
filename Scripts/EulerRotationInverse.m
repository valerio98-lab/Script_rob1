function [phi, theta, psi] = EulerRotationInverse(sequence, R, branch)
% EULERROTATIONINVERSE   Inverse Euler angles from rotation matrix
%   [phi, theta, psi] = EulerRotationInverse(seq, R, branch)
%     seq    : 3-char string specifying Euler sequence (e.g. 'xzy')
%     R      : 3x3 rotation matrix
%     branch : 'pos' or 'neg' solution for theta
% Returns angles [phi, theta, psi] in radians.

    % Validate inputs
    if ~ischar(sequence) || numel(sequence)~=3
        error('Sequence must be a 3-character string.');
    end
    if ~ismatrix(R) || ~all(size(R)==[3 3])
        error('R must be a 3x3 matrix.');
    end
    branch = lower(string(branch));
    if branch ~= "pos" && branch ~= "neg"
        error('branch must be ''pos'' or ''neg''');
    end
    cond = (branch == "pos");

    % Preallocate
    phi = NaN; theta = NaN; psi = NaN;

    switch lower(sequence)
        case 'xzy'
            % sin(theta) and |cos(theta)|
            sinth = -R(1,2);
            costh = sqrt(R(1,1)^2 + R(1,3)^2);
            % choose branch
            if cond
                theta = atan2(sinth,  costh);
            else
                theta = atan2(sinth, -costh);
            end
            if abs(cos(theta)) < eps
                error('Singular configuration: cos(theta) ~ 0');
            end
            % compute psi, phi
            psi = atan2( R(1,3)/cos(theta), R(1,1)/cos(theta) );
            phi = atan2( R(3,2)/cos(theta), R(2,2)/cos(theta) );

        % -- Other sequences follow the same pattern: define sinth & costh
        case 'xyz'
            sinth = R(1,3);
            costh = sqrt(R(1,1)^2 + R(1,2)^2);
            if cond
                theta = atan2(sinth,  costh);
            else
                theta = atan2(sinth, -costh);
            end
            if abs(cos(theta)) < eps
                error('Singular configuration: cos(theta) ~ 0');
            end
            psi = atan2(-R(1,2)/cos(theta), R(1,1)/cos(theta));
            phi = atan2(-R(2,3)/cos(theta), R(3,3)/cos(theta));

        case 'yzx'
            sinth = R(2,1);
            costh = sqrt(R(2,2)^2 + R(2,3)^2);
            if cond
                theta = atan2(sinth,  costh);
            else
                theta = atan2(sinth, -costh);
            end
            if abs(cos(theta)) < eps
                error('Singular configuration: cos(theta) ~ 0');
            end
            psi = atan2(-R(2,3)/cos(theta), R(2,2)/cos(theta));
            phi = atan2(-R(3,1)/cos(theta), R(1,1)/cos(theta));

        case 'zxy'
            sinth = R(3,2);
            costh = sqrt(R(3,1)^2 + R(3,3)^2);
            if cond
                theta = atan2(sinth,  costh);
            else
                theta = atan2(sinth, -costh);
            end
            if abs(cos(theta)) < eps
                error('Singular configuration: cos(theta) ~ 0');
            end
            psi = atan2(-R(3,1)/cos(theta), R(3,3)/cos(theta));
            phi = atan2(-R(1,2)/cos(theta), R(2,1)/cos(theta));

        % -- Similar implementation for the 6 proper/non-proper sequences
        case 'xyx'
            sinth =  sqrt(R(1,2)^2 + R(1,3)^2);
            costh = R(1,1);
            if cond
                theta = atan2(sinth,  costh);
            else
                theta = atan2(-sinth, costh);
            end
            if abs(sin(theta)) < eps
                error('Singular configuration: sin(theta) ~ 0');
            end
            psi = atan2( R(1,2)/sin(theta),  R(1,3)/sin(theta) );
            phi = atan2( R(2,1)/sin(theta), -R(3,1)/sin(theta) );

        case 'xzx'
            sinth =  sqrt(R(1,2)^2 + R(1,3)^2);
            costh = R(1,1);
            if cond
                theta = atan2(sinth,  costh);
            else
                theta = atan2(-sinth, costh);
            end
            if abs(sin(theta)) < eps
                error('Singular configuration: sin(theta) ~ 0');
            end
            psi = atan2( R(1,3)/sin(theta), -R(1,2)/sin(theta) );
            phi = atan2( R(3,1)/sin(theta),  R(2,1)/sin(theta) );

        case 'yxy'
            sinth =  sqrt(R(2,1)^2 + R(2,3)^2);
            costh = R(2,2);
            if cond
                theta = atan2(sinth,  costh);
            else
                theta = atan2(-sinth, costh);
            end
            if abs(sin(theta)) < eps
                error('Singular configuration: sin(theta) ~ 0');
            end
            psi = atan2( R(2,1)/sin(theta), -R(2,3)/sin(theta) );
            phi = atan2( R(1,2)/sin(theta),  R(3,2)/sin(theta) );

        case 'yzy'
            sinth =  sqrt(R(2,1)^2 + R(2,3)^2);
            costh = R(2,2);
            if cond
                theta = atan2(sinth, costh);
            else
                theta = atan2(-sinth,costh);
            end
            if abs(sin(theta)) < eps
                error('Singular configuration: sin(theta) ~ 0');
            end
            psi = atan2( R(2,3)/sin(theta),  R(2,1)/sin(theta) );
            phi = atan2( R(3,2)/sin(theta), -R(1,2)/sin(theta) );

        case 'zyz'
            sinth =  sqrt(R(3,1)^2 + R(3,2)^2);
            costh = R(3,3);
            if cond
                theta = atan2(sinth, costh);
            else
                theta = atan2(-sinth,costh);
            end
            if abs(sin(theta)) < eps
                error('Singular configuration: sin(theta) ~ 0');
            end
            psi = atan2( R(3,2)/sin(theta), -R(3,1)/sin(theta) );
            phi = atan2( R(2,3)/sin(theta),  R(1,3)/sin(theta) );

        case 'yxz'
            sinth = -R(2,3);
            costh = sqrt(R(2,1)^2 + R(2,2)^2);
            if cond
                theta = atan2(sinth, costh);
            else
                theta = atan2(sinth,-costh);
            end
            if abs(cos(theta)) < eps
                error('Singular configuration: cos(theta) ~ 0');
            end
            psi = atan2( R(2,1)/cos(theta),  R(2,2)/cos(theta) );
            phi = atan2( R(1,3)/cos(theta), -R(3,3)/cos(theta) );

        otherwise
            error('Unsupported sequence "%s".', sequence);
    end
end
