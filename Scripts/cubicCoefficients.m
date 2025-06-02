function [q, pathCoeffs, vel_t_max, vel_val_max, acc_t_max, acc_val_max] = cubicCoefficients(qi, qf, velocity, T, DK, pi_dot, pf_dot, N)
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
        %[q, C, tV, vMax, tA, aMax] = cubicCoefficients([-0.78;1.2; 1.5], [ 0.40; -0.4; 0], false, 2.0, []);

        % table(tV, vMax, tA, aMax)
        % disp(q); 
        
        %% Caso con DH:
        % f_q = [q2*cos(q1)+0.6*cos(q1+q3);
        %     q2*sin(q1)+0.6*sin(q1+q3);
        %     q1+q3]; 
        % 
        % [q, C, tV, vMax, tA, aMax] = cubicCoefficients([-0.78;1.2; 1.5], [ 0.40; -0.4; 0], true, 2.0, f_q, [0,0,0], [1,6,8]);
        % 
        % table(tV, vMax, tA, aMax)
        
        %% SENZA DH MA VELOCITà NON REST
        % [q, C, tV, vMax, tA, aMax] = cubicCoefficients([-0.78;1.2; 1.5], [ 0.40; -0.4; 0], true, 2.0, [], [0,0,0], [1,6,8]);
        % 
        % table(tV, vMax, tA, aMax)

%%

    if nargin < 8 || isempty(N),        N = 100;       end
    if nargin < 7,                      pf_dot = [];   end
    if nargin < 6,                      pi_dot = [];   end
    if nargin < 5,                      DK     = [];   end
    if nargin < 4 || isempty(T),        T      = 0;    end
    if nargin < 3 || isempty(velocity), velocity = false; end


    %--------------------------------------
    % 1) Check that qi, qf have the same length
    %--------------------------------------
    % 1) Check dimensional consistency on positions --------------------------
    num_joints = numel(qi);
    if numel(qf) ~= num_joints
        error('qi and qf must have the same length.');
    end
    
    % 2) Symbolic joint vector ------------------------------------------------
    joint_vars = sym('q', [1 num_joints]);
    
    % 3) Row-vectors di configurazione (servono comunque) ---------------------
    qi_row = qi(:).';          %  1 × N
    qf_row = qf(:).';
    
    % 4) Gestione DK / Jacobiano / velocità -----------------------------------
    if isempty(DK)                 % ****** CASO A: DK assente ******
        v_j_i = [];  v_j_f = [];
    
        if velocity
            if isempty(pi_dot), pi_dot = zeros(num_joints,1); end
            if isempty(pf_dot), pf_dot = zeros(num_joints,1); end
            if numel(pi_dot)~=num_joints || numel(pf_dot)~=num_joints
                error('pi_dot / pf_dot must have the same length as qi.');
            end
            q_dot_i = pi_dot(:);
            q_dot_f = pf_dot(:);
        else
            q_dot_i = zeros(num_joints,1);
            q_dot_f = zeros(num_joints,1);
        end
    
    else                            % ****** CASO B: DK presente ******
        J = jacobian(DK, joint_vars);
    
        v_j_i = double( subs(J, joint_vars, qi_row) );
        v_j_f = double( subs(J, joint_vars, qf_row) );
    
        if velocity
            if isempty(pi_dot), pi_dot = zeros(size(DK)); end
            if isempty(pf_dot), pf_dot = zeros(size(DK)); end
            if numel(pi_dot)~=size(DK,1) || numel(pf_dot)~=size(DK,1)
                error('pi_dot / pf_dot must have length equal to rows of DK.');
            end
    
            if size(J,1) == size(J,2)
                J_inv = inv(J);
            else
                J_inv = pinv(J);
            end
            J_inv_i = double( subs(J_inv, joint_vars, qi_row) );
            J_inv_f = double( subs(J_inv, joint_vars, qf_row) );
    
            q_dot_i = J_inv_i * pi_dot(:);
            q_dot_f = J_inv_f * pf_dot(:);
        else
            q_dot_i = zeros(num_joints,1);
            q_dot_f = zeros(num_joints,1);
        end
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

    vel_t_max     = zeros(num_joints,1);
    vel_val_max   = zeros(num_joints,1);
    acc_t_max     = zeros(num_joints,1);
    acc_val_max   = zeros(num_joints,1);


    for k = 1:num_joints
        A0 = pathCoeffs(k,1);
        A1 = pathCoeffs(k,2);
        A2 = pathCoeffs(k,3);
        A3 = pathCoeffs(k,4);
    
        % ---- vel -----
        cand_s = [0 1];
        if abs(A3) > eps
            s_star = -A2/(3*A3);
            if s_star>=0 && s_star<=1
                cand_s(end+1) = s_star;
            end
        end
        vel_s = A1 + 2*A2*cand_s + 3*A3*cand_s.^2;
    
        if T>0
            [vmax, idx]    = max(abs(vel_s)/T);
            vel_t_max(k)   = cand_s(idx)*T;
        else
            [vmax, idx]    = max(abs(vel_s));     % tempo normalizzato
            vel_t_max(k)   = cand_s(idx);
        end
        vel_val_max(k) = vmax;
    
        % ---- acc -----
        bordi = [0 1];                  % aggiungi subito prima di usarlo
    
        acc_s = 2*A2 + 6*A3*bordi;      % usa 'bordi' anche qui
        if T>0
            [amax, idy]  = max(abs(acc_s)/T^2);
            acc_t_max(k) = bordi(idy)*T;     % istante in secondi
        else
            [amax, idy]  = max(abs(acc_s));
            acc_t_max(k) = bordi(idy);       % istante in s normalizzato
        end
        acc_val_max(k) = amax;
    end


    % 9) (Optional) Pretty‐print each polynomial’s numeric coefficients:
    fprintf('\n# Max |velocity| and |acceleration|\n');
    for k = 1:num_joints
        fprintf('Joint %d:  |q̇|max = %.4g at t = %.4g s   |q̈|max = %.4g at t = %.4g s\n', ...
                k, vel_val_max(k), vel_t_max(k), acc_val_max(k), acc_t_max(k));
    end

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

    if ~isempty(v_j_i)
        disp('Velocità di giunto iniziale v_j_i (Jacobiana @ qi):');
        disp(v_j_i);
        disp('Velocità di giunto finale   v_j_f (Jacobiana @ qf):');
        disp(v_j_f);
    end

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
