function plotTaskJointTraj(p_i, p_f, A_max, V_max, scale, N)
% plotTaskJointTraj - Plot Cartesian and Joint space trajectories
%   plotTaskJointTraj(p_i, p_f, A_max, V_max)
%   plotTaskJointTraj(..., scale)
%   plotTaskJointTraj(..., scale, N)
%
% Inputs:
%   p_i   : [2×1] start position in task-space [x; y]
%   p_f   : [2×1] end position in task-space [x; y]
%   A_max : maximum task-space acceleration (scalar)
%   V_max : maximum task-space velocity (scalar)
%   scale : (optional) time scaling factor (>0, default=1)
%   N     : (optional) number of samples per segment (default=300)
%
% This function:
% 1) Computes a rest-to-rest trapezoidal profile along the line
%    between p_i and p_f using Rest2RestTraj.
% 2) Maps the scalar profile s(t) into Cartesian signals x, y, vx, vy, ax, ay.
% 3) Uses inverse kinematics and the geometric Jacobian to compute
%    joint-space trajectories q1, q2, q1_dot, q2_dot, q1_ddot, q2_ddot.
% 4) Plots position, velocity, and acceleration in both Cartesian and Joint spaces.

if nargin < 5 || isempty(scale)
  scale = 1;
end
if nargin < 6 || isempty(N)
  N = 300;
end
% Direction and length
delta_p = p_f - p_i;
L = norm(delta_p);
u = delta_p / L;
% Compute scalar profile
[time, s, s_dot, s_ddot] = Rest2RestTraj(L, A_max, V_max, scale);
% Map to task-space
x   = p_i(1) + u(1)*s;
y   = p_i(2) + u(2)*s;
vx  = u(1)*s_dot;
vy  = u(2)*s_dot;
ax  = u(1)*s_ddot;
ay  = u(2)*s_ddot;
% Preallocate joint trajectories
iN = numel(time);
q1       = zeros(1,iN);
q2       = zeros(1,iN);
q_dot    = zeros(2,iN);
% Compute joint trajectories
for i = 1:iN
  % inverse kinematics
  q1(i) = atan2(y(i), x(i));
  q2(i) = hypot(x(i), y(i));
  % Jacobian
  Jv = [-q2(i)*sin(q1(i)), cos(q1(i));
        q2(i)*cos(q1(i)), sin(q1(i))];
  % joint velocities
  q_dot(:,i) = Jv \ [vx(i); vy(i)];
end
q1_dot = q_dot(1,:);
q2_dot = q_dot(2,:);
% joint accelerations via numerical derivative
q1_ddot = gradient(q1_dot, time);
q2_ddot = gradient(q2_dot, time);
% Plotting
figure('Name','Trajectory: Cartesian vs Joint','Units','normalized','Position',[.1 .1 .8 .8]);
% Cartesian
subplot(3,2,1);
plot(time, x, 'b', time, y, 'r', 'LineWidth',1.5);
legend('x(t)','y(t)'); title('Cartesian Position'); xlabel('t [s]'); grid on;
subplot(3,2,3);
plot(time, vx, 'b', time, vy, 'r', 'LineWidth',1.5);
legend('vx(t)','vy(t)'); title('Cartesian Velocity'); xlabel('t [s]'); grid on;
subplot(3,2,5);
plot(time, ax, 'b', time, ay, 'r', 'LineWidth',1.5);
legend('ax(t)','ay(t)'); title('Cartesian Acceleration'); xlabel('t [s]'); grid on;
% Joint
subplot(3,2,2);
plot(time, q1, 'b', time, q2, 'r', 'LineWidth',1.5);
legend('q_1(t)','q_2(t)'); title('Joint Position'); xlabel('t [s]'); grid on;
subplot(3,2,4);
plot(time, q1_dot, 'b', time, q2_dot, 'r', 'LineWidth',1.5);
legend('q_1-dot','q_2-dot'); title('Joint Velocity'); xlabel('t [s]'); grid on;
subplot(3,2,6);
plot(time, q1_ddot, 'b', time, q2_ddot, 'r', 'LineWidth',1.5);
legend('q_1-ddot','q_2-ddot'); title('Joint Acceleration'); xlabel('t [s]'); grid on;
end
