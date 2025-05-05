function [time, position, velocity, acceleration] = Rest2RestTraj(L, A_max, V_max, scale)
% RESTTORESTMOTION – profilo rest-to-rest bang–coast–bang con scalatura
%
%   [t, s, v, a] = restToRestMotion(L, A_max, V_max, scale)
%
% Input:
%   L      : spostamento totale (m o rad)
%   A_max  : accelerazione massima |a(t)| ≤ A_max
%   V_max  : velocità massima |v(t)| ≤ V_max
%   scale  : fattore di scalatura temporale (>1 rallenta, <1 accelera)
%
% Output:x
%   time        : vettore temporale 0→T_scaled
%   position    : s(t)
%   velocity    : v(t)
%   acceleration: a(t)
%
  if nargin<4 || isempty(scale)
    scale = 1;  % default nessuna scalatura
  end

  % 1) calcolo Ts e T (originali)
  Ts0 = V_max / A_max;                  % Ts = V_max/A_max
  if L >= V_max^2 / A_max
    % profilo trapezoidale
    Ts = Ts0;
    T  = (L*A_max + V_max^2) / (A_max * V_max);  % T = (L·A_max+V_max^2)/(A_max·V_max)
    type = 'Trapezoidale';
  else
    % profilo triangolare
    Ts = sqrt(L / A_max);             % Ts = sqrt(L/A_max)
    V_max = A_max * Ts;               % V_peak
    T  = 2 * Ts;                      % T = 2·Ts
    type = 'Triangolare';
  end

  % 2) tempi scalati
  Ts_scaled = Ts * scale;
  T_scaled  = T  * scale;

  fprintf('Profilo %s (scaled ×%.2g): Ts=%.4f→%.4f  T=%.4f→%.4f  V_pk=%.4f\n', ...
          type, scale, Ts, Ts_scaled, T, T_scaled, V_max);

  % 3) griglie originali e scalate
  N = 300;
  t1o = linspace(0, Ts,     N);
  t2o = linspace(Ts, T-Ts,  N);
  t3o = linspace(T-Ts, T,   N);
  t_o = [t1o, t2o, t3o];             % tempo parametrico originale
  time = t_o * scale;                % tempo scalato

  % 4) profilo a(t) originale e scalato
  a_o = [ +A_max*ones(1,N), ...      % accel.
          zeros(1,N), ...            % coast
          -A_max*ones(1,N) ];        % decel.
  acceleration = a_o / scale^2;      % a_scaled(t) = a_o(t/scale)/scale^2

  % 5) profilo v(t) originale e scalato
  v1o = A_max .* t1o;
  v2o = V_max  * ones(1,N);
  v3o = A_max .* (T - t3o);
  v_o = [v1o, v2o, v3o];
  velocity = v_o / scale;            % v_scaled = v_o/scale

  % 6) profilo s(t) originale
  s1 = 0.5 * A_max .* (t1o.^2);
  s2 = V_max .* t2o - V_max^2/(2*A_max);
  s3 = -0.5 * A_max .* ((t3o - T).^2) + V_max*T - V_max^2/A_max;
  position = [s1, s2, s3];           % s_scaled = s_o (invariante)

  % controllo fine corsa
  if abs(position(end) - L) > 1e-6
    warning('s(T) = %.6f ≠ L = %.6f', position(end), L);
  end

  % 7) plot con colori
  figure('Name','Rest-to-Rest Motion (scaled)','Units','normalized','Position',[.2 .3 .5 .5]);
  subplot(3,1,1);
    plot(time, position, 'b', 'LineWidth',1.5);
    title('s(t) = posizione'); xlabel('t [s]'); ylabel('s(t)'); grid on;
  subplot(3,1,2);
    plot(time, velocity, 'r', 'LineWidth',1.5);
    title('v(t) = velocità');    xlabel('t [s]'); ylabel('v(t)'); grid on;
  subplot(3,1,3);
    plot(time, acceleration, 'g', 'LineWidth',1.5);
    title('a(t) = accelerazione');xlabel('t [s]'); ylabel('a(t)'); grid on;

end
