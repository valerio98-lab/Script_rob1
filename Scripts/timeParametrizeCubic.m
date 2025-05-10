function [t, q, dq, ddq, s, ds] = timeParametrizeCubic(pathCoeffs, timeCoeffs, T, N, doPlot)
% TIMEPARAMETRIZE  Compose geometric path q(s) with timing law s(t).

% [qSym, pathCoeffs] = cubicCoefficients(qi, qf, false, 2, DK);
% timeCoeffs = [0 0 3/4 -1/4];
% [t, q, dq, ddq, s, ds] = timeParametrizeCubic(pathCoeffs, timeCoeffs, 2, 100, false);


%   [t, q, dq, ddq] = timeParametrize(P, S, T)  – default 100 samples,
%   plots q(t), dq(t), ddq(t) in a single window if doPlot=true.
%
%   pathCoeffs : nJ×4 [a0 a1 a2 a3]
%   timeCoeffs : 1×4  [b0 b1 b2 b3]
%   T          : total duration (s)
%   N          : samples (opt, default 100)
%   doPlot     : logical, default true
%
%   Output arrays are nJ×N.
% -----------------------------------------------------------


if nargin < 4 || isempty(N),       N = 100; end
if nargin < 5 || isempty(doPlot), doPlot = true; end

[nJ, cols] = size(pathCoeffs);
if cols ~= 4, error('pathCoeffs must be nJ×4'); end
if numel(timeCoeffs) ~= 4, error('timeCoeffs must be 1×4'); end
validateattributes(T, {'numeric'}, {'scalar','>',0});

timeCoeffs = timeCoeffs(:).';

% ---- time base ----
 t  = linspace(0, T, N);
 b0 = timeCoeffs(1);
 b1 = timeCoeffs(2);
 b2 = timeCoeffs(3);
 b3 = timeCoeffs(4);
 s   = b0 + b1*t + b2*t.^2 + b3*t.^3;
 ds  =       b1     + 2*b2*t + 3*b3*t.^2;
 dds =                2*b2    + 6*b3*t;

q   = zeros(nJ,N);
dq  = zeros(nJ,N);
ddq = zeros(nJ,N);

for j = 1:nJ
    a0 = pathCoeffs(j,1); a1 = pathCoeffs(j,2); a2 = pathCoeffs(j,3); a3 = pathCoeffs(j,4);
    q_s   = a0 + a1*s + a2*s.^2 + a3*s.^3;
    dq_s  =       a1     + 2*a2*s + 3*a3*s.^2;
    ddq_s =                2*a2     + 6*a3*s;
    q(j,:)   = q_s;
    dq(j,:)  = dq_s .* ds;
    ddq(j,:) = ddq_s .* ds.^2 + dq_s .* dds;
end

% ---- plotting ----
if doPlot
    figure('Name','Joint trajectories');
    tl = tiledlayout(3,1, 'TileSpacing','Compact');
    nexttile; plot(t,q,'LineWidth',1.2);  ylabel('q');  grid on; legend(compose('q_%d',1:nJ));
    nexttile; plot(t,dq,'LineWidth',1.2); ylabel('dq'); grid on;
    nexttile; plot(t,ddq,'LineWidth',1.2); xlabel('Time [s]'); ylabel('ddq'); grid on;
    title(tl,'Joint Position, Velocity, Acceleration vs Time');
end