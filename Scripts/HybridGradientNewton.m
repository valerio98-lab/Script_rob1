function [q_out, info] = HybridGradientNewton(q_sym, desired_point, direct_kyn, q0, varargin)
%HYBRIDGRADIENTNEWTON  Inverse kinematics: gradient warm‑start then Newton refinement
% 
% f_q = [0.5*cos(q1) + 0.5*cos(q1+q2)*cos(q3);
%     0.5*sin(q1)+0.5*sin(q1+q2)*cos(q3);
%     0.5+0.5*sin(q3)];
% 
% 
% [q_out, info] = HybridGradientNewton([q1 q2 q3], [0.3;-0.3;0.7], f_q, [-pi/4 pi/4 pi/4], ...
%     'CartesianNewtonTol', 10e-4, 'Verbose', true);
% 
% disp(info.newton)
% disp(info.grad)
%
%   [q_out, info] = HybridGradientNewton(q_sym, desired_point, direct_kyn, q0, ...)
%
% Phase‑1  GradientMethod  (coarse tolerance / limited iterations)
% Phase‑2  NewtonMethod    (fine tolerance)
% If Phase‑1 does not meet the switching tolerance, Phase‑2 can still be
% attempted if you set 'AlwaysNewton', true (default).
%
% INPUTS (positional)
%   q_sym          1×n symbolic vector of joints
%   desired_point  m×1 numeric target pose
%   direct_kyn     m×1 symbolic f(q)
%   q0             1×n numeric initial guess
%
% OPTIONAL NAME–VALUE PAIRS
%   'GradTol'       1e-3    – switch to Newton when ||e||<GradTol
%   'GradMaxIter'   300     – max iterations for GradientMethod
%   'CartesianNewtonTol'     1e-8    – CartesianTolerance for NewtonMethod
%   'AlwaysNewton'  true    – still run Newton even if GradTol not reached
%   'GradOpts'      struct  – extra args to GradientMethod
%   'NewtonOpts'    struct  – extra args to NewtonMethod
%   'Verbose'       false   – high‑level log
%
% OUTPUTS
%   q_out                n×1 solution (Newton if successful, else best grad)
%   info                 struct
%       .grad            – GradientMethod log (see GradientMethod)
%       .newton          – NewtonMethod log  (empty if not run)
%       .phase1Converged – bool
%       .phase2Converged – bool (true only if Newton ran & met tol)
%       .totalIter       – total iterations (grad + newton)
%
% DEPENDS ON GradientMethod.m  NewtonMethod.m
% --------------------------------------------------------------------
% Matteo, 2025 – revised after user feedback

% ------- parse wrapper options ---------------------------------------
p = inputParser;
addParameter(p,'GradTol',      1e-3,  @(x) isnumeric(x)&&isscalar(x));
addParameter(p,'GradMaxIter',  300,   @(x) isnumeric(x)&&isscalar(x));
addParameter(p,'CartesianNewtonTol',    1e-8,  @(x) isnumeric(x)&&isscalar(x));
addParameter(p,'AlwaysNewton', true,  @(x) islogical(x)&&isscalar(x));
addParameter(p,'GradOpts',     struct, @(x) isstruct(x));
addParameter(p,'NewtonOpts',   struct, @(x) isstruct(x));
addParameter(p,'Verbose',      false, @(x) islogical(x)&&isscalar(x));
parse(p,varargin{:});
W = p.Results;

% ---------- Phase 1 : Gradient descent --------------------------------
GradOpts        = W.GradOpts;
GradOpts.MaxIterations      = W.GradMaxIter;
GradOpts.CartesianTolerance = W.GradTol;

[q_grad, guesses_g, fk_hist_g, err_g, alpha_hist] = ...
        GradientMethod(q_sym, desired_point, direct_kyn, q0, GradOpts);

phase1Conv = err_g(end) < W.GradTol;
if W.Verbose
    fprintf('[Hybrid] Phase‑1 (Gradient) finished: err = %.3e after %d iter (target %.1e)\n', ...
            err_g(end), numel(err_g), W.GradTol);
end

% Log grad phase
info.grad.q_out         = q_grad;
info.grad.guesses       = guesses_g;
info.grad.fk_history    = fk_hist_g;
info.grad.cartesian_err = err_g;
info.grad.alpha_hist    = alpha_hist;
info.phase1Converged    = phase1Conv;

% ---------- Decide whether to run Newton ------------------------------
runNewton = phase1Conv || W.AlwaysNewton;

info.newton = [];
phase2Conv  = false;
q_refine    = q_grad;

if runNewton
    NewtonOpts                         = W.NewtonOpts;
    NewtonOpts.CartesianTolerance      = W.CartesianNewtonTol;
    NewtonOpts.MaxIterations           = getfield_with_default(NewtonOpts,'MaxIterations',1000);

    try
        [q_new, guesses_n, fk_hist_n, err_n] = ...
            NewtonMethod(q_sym, desired_point, direct_kyn, q_refine, NewtonOpts);

        info.newton.q_out         = q_new;
        info.newton.guesses       = guesses_n;
        info.newton.fk_history    = fk_hist_n;
        info.newton.cartesian_err = err_n;

        phase2Conv = err_n(end) < W.CartesianNewtonTol;

        if W.Verbose
            fprintf('[Hybrid] Phase‑2 (Newton) finished: err = %.3e after %d iter (target %.1e)\n', ...
                    err_n(end), numel(err_n), W.CartesianNewtonTol);
        end
    catch ME
        warning('[Hybrid] NewtonMethod error: %s', ME.message);
    end
else
    if W.Verbose
        fprintf('[Hybrid] Newton phase skipped.\n');
    end
end

info.phase2Converged = phase2Conv;

% ---------- choose output q ------------------------------------------
if phase2Conv
    q_out = info.newton.q_out;
else
    q_out = q_grad;   % best available
end

% ---------- total iteration count ------------------------------------
info.totalIter = numel(err_g);
if ~isempty(info.newton)
    info.totalIter = info.totalIter + size(info.newton.guesses,1);
end
end

% -------- helper -----------------------------------------------------
function v = getfield_with_default(s,f,def)
    if isfield(s,f); v = s.(f); else; v = def; end
end
