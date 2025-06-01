function [q, pathCoeffs] = cubicCoefficients(qi, qf, velocity, T, DK, pi_dot, pf_dot, N)
    % CUBICCOEFFICIENTS  Joint cubic interpolation *and* Cartesian path plot.
    %
    %   q = cubicCoefficients(qi,qf,velocity,T,DK,pi_dot,pf_dot, N) uses the
    %   symbolic direct-kinematics vector DK(q) to build third-order polynomials
    %   that move the manipulator from configuration qi to qf in normalized time
    %   s ∈ [0,1] (physical time t = s·T). If T>0 it plots:
    %       • joint trajectories (as before)
    %       • Cartesian path of the end-effector in 3-D (or 2-D if DK is 2×1)
    %
    %   N is the number of points for plotting (default = 100).
       %% REST 2 REST
%     [q, coeffs] = cubicCoefficients(-.78, 0, false, 1,[]);
%       disp(q); 
%%

    if nargin < 8 || isempty(N),        N = 100;       end
    if nargin < 7,                      pf_dot = [];   end
    if nargin < 6,                      pi_dot = [];   end
    if nargin < 5,                      DK     = [];   end
    if nargin < 4 || isempty(T),        T      = 0;    end
    if nargin < 3 || isempty(velocity), velocity = false; end


    %--------------------------------------
    % 1) Check that qi, qf have the same length
    num_joints_i = length(qi);
    num_joints_f = length(qf);
    if num_joints_i ~= num_joints_f
        error('The number of initial joints does not match the number of final joints.');
    end
    num_joints = num_joints_i;

    % 2) Build the symbolic joint‐variable vector as a 1×N row:
    joint_vars = sym('q', [1 num_joints]);   % e.g. [q1, q2, ..., qN]

    % 3) Compute the Jacobian symbolically:
    J = jacobian(DK, joint_vars);

    %-----------------------------
    % 4) Convert qi, qf into 1×N row vectors BEFORE calling subs:
    qi_row = qi(:).';   % now size = [1×N]
    qf_row = qf(:).';   % now size = [1×N]

    % 5) Substitute into J, then convert to double for numeric use
    v_j_i = double(subs(J, joint_vars, qi_row));   % J evaluated at qi
    v_j_f = double(subs(J, joint_vars, qf_row));   % J evaluated at qf

    if velocity
        % If DK is square, use inv; otherwise pinv:
        if isempty(DK)                       % <<< NUOVO CASO >>>
            %  ➜  pi_dot e pf_dot sono GIA' velocità di giunto
            q_dot_i = pi_dot(:);
            q_dot_f = pf_dot(:);
        else

            if size(J,1) == size(J,2)
                J_inv = inv(J);
            else
                J_inv = pinv(J);
            end
            % Evaluate J_inv at qi_row and qf_row, then convert to double:
            J_inv_i = double(subs(J_inv, joint_vars, qi_row));
            J_inv_f = double(subs(J_inv, joint_vars, qf_row));
    
            % Multiply by Cartesian velocity vectors pi_dot, pf_dot:
            q_dot_i = J_inv_i * pi_dot(:);
            q_dot_f = J_inv_f * pf_dot(:);
        end
    else
        q_dot_i = zeros(num_joints,1);
        q_dot_f = zeros(num_joints,1);
    end

    %-----------------------------------------
    % 6) Now reshape qi, qf back into N×1 columns for the polynomial math:
    qi = qi(:);
    qf = qf(:);

    % Check consistency
    if ~( all(size(qi) == size(qf)) && ...
          all(size(qi) == size(q_dot_i)) && ...
          all(size(qi) == size(q_dot_f)) )
        error('Dimension mismatch among qi, qf, q_dot_i, q_dot_f');
    end

    % (In your original you force specific q_dot values—for testing—here’s that snippet:)
    % q_dot_i = [-2.5; 2.5];
    % q_dot_f = [-0.3; -0.1];

    %-----------------------------------------
    % 7) Build cubic‐polynomial coefficients:
    delta_q = qf - qi;
    a0      = qi;
    a1      = q_dot_i;
    a2      = 3*delta_q - (q_dot_f + 2*q_dot_i);
    a3      = -2*delta_q + (q_dot_f + q_dot_i);

    % pathCoeffs is an N×4 matrix, where each row i = [a0(i), a1(i), a2(i), a3(i)]
    pathCoeffs = [a0, a1, a2, a3];

    %-----------------------------------------
    % 8) Return symbolic expressions for q_i(s):
    syms s real
    q = cell(1, num_joints);
    for k = 1:num_joints
        q{k} = a0(k) + a1(k)*s + a2(k)*s^2 + a3(k)*s^3;
    end

    % 9) (Optional) Pretty‐print each polynomial’s numeric coefficients:
    for k = 1:num_joints
        pc = pathCoeffs(k,:);
        txt = sprintf('Poly %d:',k);
        for d = 0:3                      % gradi 0→3
            if abs(pc(d+1)) > 1e-12      % soglia piccola
                txt = sprintf('%s  %+g*s^%d',txt,pc(d+1),d);
            end
        end
        disp(txt)
    end

    disp('Velocità di giunto iniziale v_j_i (Jacobiana @ qi):');
    disp(v_j_i);
    disp('Velocità di giunto finale   v_j_f (Jacobiana @ qf):');
    disp(v_j_f);

    %-----------------------------------------
    % 10) If T > 0, plot joint trajectories and Cartesian path:
    if T > 0
        plotTrajectories(a0, a1, a2, a3, T);
        if ~isempty(DK)
            sVals = linspace(0,1,N);
            qTraj = a0 + a1.*sVals + a2.*(sVals.^2) + a3.*(sVals.^3);
            qTraj = reshape(qTraj, num_joints, N);

            % Evaluate the Cartesian path point‐by‐point:
            cart = zeros(length(DK), N);
            for idx = 1:N
                cart(:,idx) = double(subs(DK, joint_vars, qTraj(:,idx).'));
            end

            figure('Name','Cartesian path');
            if size(cart,1) >= 3
                plot3(cart(1,:), cart(2,:), cart(3,:), 'LineWidth', 2);
                xlabel('X'); ylabel('Y'); zlabel('Z'); grid on; axis equal
                title('End-effector Cartesian Path');
                setAxisWithPad(cart(1,:), cart(2,:), cart(3,:), 0.1);
            else
                plot(cart(1,:), cart(2,:), 'LineWidth', 2);
                xlabel('X'); ylabel('Y'); grid on; axis equal
                title('End-effector Cartesian Path');
                setAxisWithPad(cart(1,:), cart(2,:), [], 0.20);
            end
        else
            warning('DK is empty – cannot compute Cartesian path.');
        end
    end
end  % ---- end cubicCoefficients function ----


% =====================================================================
function setAxisWithPad(x, y, z, pad)
    % If 2-D:
    if isempty(z)
        xr = max(x) - min(x);
        yr = max(y) - min(y);
        if xr == 0, xr = 1; end
        if yr == 0, yr = 1; end
        xlim([min(x) - pad*xr, max(x) + pad*xr]);
        ylim([min(y) - pad*yr, max(y) + pad*yr]);
    else
        % 3-D case:
        xr = max(x) - min(x);
        yr = max(y) - min(y);
        zr = max(z) - min(z);
        if xr == 0, xr = 1; end
        if yr == 0, yr = 1; end
        if zr == 0, zr = 1; end
        xlim([min(x) - pad*xr, max(x) + pad*xr]);
        ylim([min(y) - pad*yr, max(y) + pad*yr]);
        zlim([min(z) - pad*zr, max(z) + pad*zr]);
    end
end
