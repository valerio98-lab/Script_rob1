function [q, pathCoeffs] = cubicCoefficients(qi, qf, velocity, T, DK, pi_dot, pf_dot, N)
% CUBICCOEFFICIENTS  Joint cubic interpolation *and* Cartesian path plot.
%   q = cubicCoefficients(qi,qf,velocity,T,pi_dot,pf_dot,DK) uses the
%   symbolic direct‑kinematics vector DK(q) to build third‑order polynomials
%   that move the manipulator from configuration qi to qf in normalised time
%   s∈[0,1] (physical time t = s·T). If T>0 it plots:
%       • joint trajectories (as before)
%       • Cartesian path of the end‑effector in 3‑D (or 2‑D if DK is 2×1)
%
%   q = cubicCoefficients(..., N) lets you choose the sample resolution for
%   the plots (default N = 100).
%


if nargin < 8 || isempty(N);  N = 100;  end
if nargin<4 || isempty(T);    T = 0;   end
if nargin < 4 && (velocity || isempty(velocity))
    error('Check the input parameters');
end

num_joints_i = length(qi);
num_joints_f = length(qf);

if num_joints_i ~= num_joints_f
    error('The number of Initial joints does not match final ones');
end
joint_vars = sym('q', [1 num_joints_i]);
J = jacobian(DK, joint_vars);

%calcolo le velocità di giunto iniziali e finali 

v_j_i = subs(J, joint_vars, qi);
v_j_f = subs(J, joint_vars, qf);

if velocity
    if size(J,1) == size(J,2)
        J_inv = inv(J);
    else
        J_inv = pinv(J);
    end
    J_inv_i = double(subs(J_inv, joint_vars, qi));
    J_inv_f = double(subs(J_inv, joint_vars, qf));
    q_dot_i = J_inv_i * pi_dot(:);
    q_dot_f = J_inv_f * pf_dot(:);
else
    q_dot_i = zeros(num_joints_i,1);
    q_dot_f = zeros(num_joints_i,1);
end

qi      = qi(:);
qf      = qf(:);

if ~( all(size(qi)==size(qf)) && all(size(qi)==size(q_dot_i)) && all(size(qi)==size(q_dot_f)) )
    error('Dimension mismatch among qi, qf, q_dot_i, q_dot_f');
end

q_dot_i = [-2.5;2.5];
q_dot_f = [-0.3;-0.1];

delta_q = qf - qi;
a0      = qi;
a1      = q_dot_i;
a2      = 3*delta_q - (q_dot_f + 2*q_dot_i);
a3      = -2*delta_q + (q_dot_f + q_dot_i);

pathCoeffs = [a0 a1 a2 a3]; %nJ x 4 numeric 


syms s real
% qSym = arrayfun(@(k) a0(k) + a1(k)*s + a2(k)*s^2 + a3(k)*s^3, 1:num_joints_i, 'uni',0);
q = cell(1, num_joints_i);
for i = 1:num_joints_i
    q{i} = a0(i) + a1(i)*s + a2(i)*s^2 + a3(i)*s^3;
end

% --- pretty print coefficients ---
for i = 1:numel(q)
    expr = expand(q{i});
    c_desc = sym2poly(expr);
    fprintf('Poly %d: %.6f + %.6f*s + %.6f*s^2 + %.6f*s^3\n', i, c_desc(end), c_desc(end-1), c_desc(end-2), c_desc(end-3));
end


disp('Velocità di giunto iniziale v0 con Jacobiana valuta in qi');
disp(v_j_i)
disp('Velocità di giunto iniziale vf con Jacobiana valuta in qf');
disp(v_j_f)

% -------- plotting --------
if T > 0
    plotTrajectories(a0,a1,a2,a3,T);
    if ~isempty(DK)
        sVals = linspace(0,1,N);
        qTraj = a0 + a1.*sVals + a2.*(sVals.^2) + a3.*(sVals.^3);
        qTraj = reshape(qTraj, num_joints_i, N);
        cart  = zeros(length(DK), N);
        for k = 1:N
            cart(:,k) = double(subs(DK, sym('q',[1 num_joints_i]), qTraj(:,k)'));
        end
        figure('Name','Cartesian path');
        if size(cart,1) >= 3
            plot3(cart(1,:), cart(2,:), cart(3,:), 'LineWidth', 2);
            xlabel('X'); ylabel('Y'); zlabel('Z'); grid on; axis equal
            title('End‑effector Cartesian Path');
            % ------ fixed padded limits ------
            setAxisWithPad(cart(1,:), cart(2,:), cart(3,:), PAD);
        else
            plot(cart(1,:), cart(2,:), 'LineWidth', 2);
            xlabel('X'); ylabel('Y'); grid on; axis equal
            title('End‑effector Cartesian Path');
            setAxisWithPad(cart(1,:), cart(2,:), [], 0.20);
        end
    else
        warning('DK is empty – cannot compute Cartesian path.');
    end
end
end  % ---- main function ----

% =====================================================================
function setAxisWithPad(x, y, z, pad)
%SETAXISWITHPAD  Apply fixed relative margin PAD to axis limits.
    if isempty(z)  % 2‑D case
        xr = max(x)-min(x); yr = max(y)-min(y);
        if xr==0, xr=1; end; if yr==0, yr=1; end
        xlim([min(x)-pad*xr, max(x)+pad*xr]);
        ylim([min(y)-pad*yr, max(y)+pad*yr]);
    else           % 3‑D case
        xr=max(x)-min(x); yr=max(y)-min(y); zr=max(z)-min(z);
        if xr==0, xr=1; end; if yr==0, yr=1; end; if zr==0, zr=1; end
        xlim([min(x)-pad*xr, max(x)+pad*xr]);
        ylim([min(y)-pad*yr, max(y)+pad*yr]);
        zlim([min(z)-pad*zr, max(z)+pad*zr]);
    end
end
