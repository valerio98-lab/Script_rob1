function plot_motion(path_fun, s_range, varargin)
%PLOT_MOTION  Generic visualiser for planar robot paths.
%
%   PLOT_MOTION(PATH_FUN, S_RANGE) samples the geometric path returned by the
%   handle PATH_FUN over the parameter domain S_RANGE = [s0 s1] and produces:
%       • Cartesian plot x(s) vs y(s) (geometry)
%       • Orientation φ(s) (rad)  -- normal to the left of the tangent
%       • Curvature κ(s) (1/length)
%
%   Optional name‑value pairs let you add temporal information:
%       's2t', ST_FUN      handle mapping t → s(t)
%       't_range', [t0 t1] time domain (required when 's2t' supplied)
%
%   Additional options:
%       'orientation', ORIENT_FUN   handle φ(s) if you already have it
%       'N', N_SAMPLES              number of samples (default 500)
%
%   The function auto‑detects whether temporal plots must be generated.  If
%   ST_FUN is provided, it will also show φ(t), speed v(t), and curvature
%   κ(t).
%
%   Example
%   -------
%   R  = 0.25;        Pc = [0 0];
%   path = @(s) Pc + R*[cos(s); sin(s)];   % circle
%   st   = @(t) 2*pi*t/5;                  % 5‑s lap, uniform speed
%   plot_motion(path,[0 2*pi],'s2t',st,'t_range',[0 5]);
%


%% --- Parse input -------------------------------------------------------
p = inputParser;
p.addRequired ( 'path_fun', @(f) isa(f,'function_handle') );
p.addRequired ( 's_range',  @(v) isnumeric(v) && numel(v)==2 );
p.addParameter( 's2t',        [], @(f) isempty(f)||isa(f,'function_handle') );
p.addParameter( 't_range',    [], @(v) isempty(v)|| (isnumeric(v)&&numel(v)==2) );
p.addParameter( 'orientation',[], @(f) isempty(f)||isa(f,'function_handle') );
p.addParameter( 'N',          500, @(n) isnumeric(n)&&isscalar(n)&&n>=5 );

p.parse(path_fun,s_range,varargin{:});
opts = p.Results;

s0 = s_range(1);  s1 = s_range(2);
Ns = opts.N;

%% --- Sample path in s ---------------------------------------------------
S  = linspace(s0,s1,Ns);
P  = path_fun(S);                   % expects 2×Ns
if size(P,1)~=2
    error('PATH_FUN must return a 2×N array for planar paths.');
end
x = P(1,:);  y = P(2,:);

%% --- First and second derivatives w.r.t s (numeric) --------------------
ds    = S(2)-S(1);
xd_s  = gradient(x,ds);
yd_s  = gradient(y,ds);
xdd_s = gradient(xd_s,ds);
ydd_s = gradient(yd_s,ds);

%% --- Orientation --------------------------------------------------------
if isempty(opts.orientation)
    % tangent vector is (xd_s, yd_s) -> normal (to the left) is +pi/2
    phi_s = atan2(yd_s, xd_s) + pi/2;
else
    phi_s = opts.orientation(S);
end

%% --- Curvature κ(s) -----------------------------------------------------
kappa_s = (xd_s.*ydd_s - yd_s.*xdd_s) ./ (xd_s.^2 + yd_s.^2).^(3/2);

%% --- Plot geometry & s‑domain quantities --------------------------------
figure('Name','Path visualiser');
subplot(2,2,1);  plot(x,y,'LineWidth',1.5); axis equal; grid on;
xlabel('x'); ylabel('y'); title('Geometry');

subplot(2,2,3); plot(S,phi_s); grid on;
xlabel('s'); ylabel('\phi(s) [rad]'); title('Orientation vs s');

subplot(2,2,4); plot(S,kappa_s); grid on;
xlabel('s'); ylabel('κ(s) [1/len]'); title('Curvature vs s');

%% --- If temporal law supplied ------------------------------------------
if ~isempty(opts.s2t)
    if isempty(opts.t_range)
        error('Provide ''t_range'' when using ''s2t''.');
    end
    t0=opts.t_range(1); t1=opts.t_range(2);
    Nt=Ns;  T=linspace(t0,t1,Nt);
    S_t = opts.s2t(T);
    P_t = path_fun(S_t);
    x_t = P_t(1,:); y_t = P_t(2,:);

    % derivatives via gradient in time
    dt = T(2)-T(1);
    xd_t = gradient(x_t,dt); yd_t = gradient(y_t,dt);
    v_t  = sqrt(xd_t.^2 + yd_t.^2);
    phi_t = interp1(S,phi_s,S_t,'linear','extrap');
    kappa_t = interp1(S,kappa_s,S_t,'linear','extrap');

    subplot(2,2,2); plot(T,v_t); grid on;
    xlabel('t'); ylabel('v(t) [units/s]'); title('Speed vs time');

    figure('Name','Temporal');
    subplot(3,1,1); plot(T,phi_t); grid on;
    xlabel('t'); ylabel('\phi(t)'); title('Orientation vs time');
    subplot(3,1,2); plot(T,kappa_t); grid on;
    xlabel('t'); ylabel('κ(t)'); title('Curvature vs time');
    subplot(3,1,3); plot(x_t,y_t); axis equal; grid on;
    xlabel('x'); ylabel('y'); title('Geometry coloured by time');
end

end
