function J_inv = JacobianInversion(J)
%JACOBIANINVERSION Robust Jacobian inversion with Damped Least Squares
%   J_inv = JacobianInversion(J) returns:
%     – inv(J)                   if J è quadrata e ben condizionata
%     – pinv(J) con damping      in tutti gli altri casi, per evitare esplosioni
%
%   Parametri interni:
tol_det    = 1e-4;    % soglia minima per det(J) o det(J'*J)
tol_cond   = 1e3;     % soglia massima per il condition number
lambda     = 1e-2;    % damping factor
J = double(J);
[m,n] = size(J);
r     = rank(J);

if m==n
    % --- caso quadrato ---
    d = det(J);

    if abs(d)>tol_det && cond(J)<tol_cond
        % inversione diretta
        J_inv = inv(J);
    else
        warning('Jacobian near singular (det=%.3g, cond=%.3g). Using DLS.', d, cond(J));
        J_inv = (J'*J + lambda^2*eye(n)) \ J';
    end

elseif m>n && r==n
    % --- caso tall & full-column-rank (underdetermined) ---
    JTJ = J'*J;
    c   = cond(JTJ);
    if abs(det(JTJ))>tol_det && c<tol_cond
        J_inv = JTJ \ J';          % pseudo-inversa standard
    else
        warning('JTJ near singular (det=%.3g, cond=%.3g). Using DLS.', det(JTJ), c);
        J_inv = (JTJ + lambda^2*eye(n)) \ J';
    end

elseif m<n && r==m
    % --- caso wide & full-row-rank (overdetermined) ---
    JJT = J*J';
    c   = cond(JJT);
    if abs(det(JJT))>tol_det && c<tol_cond
        J_inv = J' / JJT;          % pseudo-inversa standard
    else
        warning('JJT near singular (det=%.3g, cond=%.3g). Using DLS.', det(JJT), c);
        J_inv = J' / (JJT + lambda^2*eye(m));
    end

else
    % --- caso rank-deficient ---
    warning('Jacobian rank-deficient (rank=%d). Using DLS pseudo-inverse.', r);
    J_inv = J' / (J*J' + lambda^2*eye(m));
end
end
